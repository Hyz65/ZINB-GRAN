from __future__ import division
from __future__ import print_function
import numpy as np
import time
import tensorflow.compat.v1 as tf
from util.utils import *
from ZINB_GRAN.models import ZINB_GRAN,Discriminator




def sample_zinb(shape, pi, r, p):
    """
    采样服从零膨胀负二项分布（ZINB）的张量。
    
    参数:
        shape: 输出张量的形状
        pi: 零膨胀概率（0 到 1 之间）
        r: 负二项分布的失败次数参数（可以是标量或与 shape 兼容的张量）
        p: 成功概率参数（0 到 1 之间，可广播）

    返回:
        一个 ZINB 分布的样本张量。
    """
    # Step 1: 是否为膨胀的零
    is_zero = tf.less(tf.random.uniform(shape), pi)
    
    # Step 2: 采样负二项分布（Negative Binomial）
    # tf.random.gamma(shape, alpha) 生成的是 gamma 分布
    # 使用 gamma-poisson 混合来采样负二项分布：
    gamma_sample = tf.random.gamma(shape, r, dtype=tf.float32) / (1 - p)
    negbin_sample = tf.random.poisson(shape=[], lam=gamma_sample, dtype=tf.float32)
    
    # Step 3: 合成 ZINB
    zinb_sample = tf.where(is_zero, tf.zeros_like(negbin_sample), negbin_sample)
    
    return zinb_sample

# 示例使用：
    #shape = tf.shape(latent)
    #zinb_sample = sample_zinb(shape, pi=0.3, r=2.0, p=0.5)




def train(FLAGS, adj, features, train_arr, test_arr, labels, AM, gene_names, TF, result_path):
    # Load data
    adj, size_gene, logits_train, logits_test, train_mask, test_mask, labels = load_data(
        adj, train_arr, test_arr, labels, AM)

    if FLAGS.model == 'ZINB_GRAN':
        model_func = ZINB_GRAN
    else:
        raise ValueError('Invalid argument for model: ' + str(FLAGS.model))

    # Disable eager mode for TF1 compatibility
    tf.compat.v1.disable_eager_execution()

    # Define placeholders
    placeholders = {
        'adjacency_matrix': tf.placeholder(tf.int32, shape=adj.shape),
        'features': tf.placeholder(tf.float32, shape=features.shape),
        'labels': tf.placeholder(tf.float32, shape=(None, logits_train.shape[1])),
        'labels_mask': tf.placeholder(tf.int32),
        'negative_mask': tf.placeholder(tf.int32)
    }

    input_dim = features.shape[1]
    model = model_func(placeholders, input_dim, size_gene, FLAGS.dim)

    # Initialize Discriminator
    discriminator = Discriminator(
        input_dim=FLAGS.dim,
        hidden_dim1=FLAGS.hidden1,
        hidden_dim2=FLAGS.dim,
        output_dim=1
    )

    # Generator latent output
    latent = model.layers[1](model.layers[0](placeholders['features']))
    #real_latent = tf.random.normal(tf.shape(latent))
    shape = tf.shape(latent)
    real_latent = sample_zinb(shape, pi=0.6, r=2.0, p=0.4)
    # Discriminator outputs
    d_real = discriminator(real_latent)
    d_fake = discriminator(latent)

    # Losses
    d_loss = -tf.reduce_mean(tf.math.log(d_real + 1e-10) + tf.math.log(1 - d_fake + 1e-10))
    g_adv_loss = -tf.reduce_mean(tf.math.log(d_fake + 1e-10))
    model._loss()
    total_g_loss = 0.1*g_adv_loss + model.loss

    # Optimizers
    d_optimizer = tf.train.AdamOptimizer(learning_rate=FLAGS.learning_rate)
    g_optimizer = tf.train.AdamOptimizer(learning_rate=FLAGS.learning_rate)

    d_train_op = d_optimizer.minimize(d_loss, var_list=discriminator.trainable_variables)
    g_train_op = g_optimizer.minimize(total_g_loss, var_list=model.layers[0].trainable_variables + model.layers[1].trainable_variables)

    # Init session
    sess = tf.Session()
    sess.run(tf.global_variables_initializer())

    def evaluate(adj, features, labels, mask, placeholders):
        t_test = time.time()
        #feed_dict_val = construct_feed_dict(adj, features, labels, mask, generate_mask(labels, FLAGS.ratio, len(test_arr), size_gene)[0], placeholders)
        feed_dict_val = construct_feed_dict(adj, features, labels, mask, negative_mask, placeholders)
        loss_val, acc_val = sess.run([model.loss, model.accuracy], feed_dict=feed_dict_val)
        return loss_val, 1 - acc_val, time.time() - t_test

    for epoch in range(FLAGS.epochs):
        t = time.time()
        negative_mask, _ = generate_mask(labels, FLAGS.ratio, len(train_arr), size_gene)
        feed_dict = construct_feed_dict(adj, features, logits_train, train_mask, negative_mask, placeholders)

        # Update Discriminator
        sess.run(d_train_op, feed_dict=feed_dict)

        # Update Generator (ZINB_GRAN + adversarial loss)
        loss_val, acc_val, _ = sess.run([model.loss, model.accuracy, g_train_op], feed_dict=feed_dict)

        print("Epoch:", '%04d' % (epoch + 1),
              "train_loss=", "{:.5f}".format(loss_val),
              "train_acc=", "{:.5f}".format(1 - acc_val),
              "time=", "{:.5f}".format(time.time() - t))

    print("Optimization Finished!")

    # Testing
    test_negative_mask, _ = generate_mask(labels, FLAGS.ratio, len(test_arr), size_gene)
    test_cost, test_acc, test_duration = evaluate(adj, features, logits_test, test_mask, placeholders)

    print("Test set results:", "cost=", "{:.5f}".format(test_cost),
          "accuracy=", "{:.5f}".format(test_acc), "time=", "{:.5f}".format(test_duration))

    feed_dict_val = construct_feed_dict(adj, features, logits_test, test_mask, test_negative_mask, placeholders)
    outs = sess.run(model.outputs, feed_dict=feed_dict_val)
    outs = np.array(outs)[:, 0].reshape((size_gene, size_gene))

    # Save predicted matrix
    logits_train = logits_train.reshape(outs.shape)
    TF_mask = np.zeros(outs.shape)
    for i, item in enumerate(gene_names):
        for j in range(len(gene_names)):
            if i == j or (logits_train[i, j] == 1):
                continue
            if item in TF:
                TF_mask[i, j] = 1

    idx_rec, idx_send = np.where(TF_mask)
    results = pd.DataFrame({
        'Gene1': np.array(gene_names)[idx_rec],
        'Gene2': np.array(gene_names)[idx_send],
        'EdgeWeight': outs[idx_rec, idx_send]
    })
    results = results.sort_values(by='EdgeWeight', ascending=False)
    return results



'''
1.从基因表达矩阵和标签数据中获得labels, adj_sparse, AM_sparse, var_names, TF, node_feat
labels 是一个包含基因对之间关系的列表，每个元素是一个三元组
2.根据labels数据将数据分成训练集和测试集
3.以adj, train_arr, test_arr, labels, AM作为关键词通过load_data获得adj, size_gene, logits_train, logits_test, train_mask, test_mask, labels
其中logits_train表示基因与基因之间的关系，类似于邻接矩阵，train_mask表示那些基因在训练中不可用 
4.随机产生一些负样本
5.将adj, features, logits_train, train_mask,negative_mask作为关键词输入模型中进行训练
'''