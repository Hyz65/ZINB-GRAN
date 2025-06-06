from ZINB_GRAN.layers import *
from ZINB_GRAN.metrics import *

flags = tf.compat.v1.flags
FLAGS = flags.FLAGS


class Model(object):
    def __init__(self, **kwargs):
        allowed_kwargs = {'name', 'logging'}
        for kwarg in kwargs.keys():
            assert kwarg in allowed_kwargs, 'Invalid keyword argument: ' + kwarg
        name = kwargs.get('name')
        if not name:
            name = self.__class__.__name__.lower()
        self.name = name

        logging = kwargs.get('logging', False)
        self.logging = logging

        self.vars = {}
        self.placeholders = {}

        self.layers = []
        self.activations = []

        self.inputs = None
        self.outputs = None
        self.hid = None

        self.loss = 0
        self.accuracy = 0
        self.optimizer = None
        self.opt_op = None
        self.prediction_labs = None

    def _build(self):
        raise NotImplementedError

    def build(self):
        """ Wrapper for _build() """
        with tf.variable_scope(self.name):
            self._build()

        # activations

        self.activations.append(self.inputs)
        for layer in self.layers:
            hidden = layer(self.activations[-1])
            self.activations.append(hidden)

        self.outputs = self.activations[-1]
        self.hid = self.activations[-2]

        # Store model variables for easy access
        variables = tf.get_collection(tf.GraphKeys.GLOBAL_VARIABLES, scope=self.name)
        self.vars = {var.name: var for var in variables}

        # Build metrics
        self._loss()
        self._accuracy()

        self.opt_op = self.optimizer.minimize(self.loss)

    def predict(self):
        pass

    def hidd(self):
        pass

    def _loss(self):
        raise NotImplementedError

    def _accuracy(self):
        raise NotImplementedError

    def save(self, sess=None):
        if not sess:
            raise AttributeError("TensorFlow session not provided.")
        saver = tf.train.Saver(self.vars)
        save_path = saver.save(sess, "tmp/%s.ckpt" % self.name)
        print("Model saved in file: %s" % save_path)

    def load(self, sess=None):
        if not sess:
            raise AttributeError("TensorFlow session not provided.")
        saver = tf.train.Saver(self.vars)
        save_path = "tmp/%s.ckpt" % self.name
        saver.restore(sess, save_path)
        print("Model restored from file: %s" % save_path)



class ZINB_GRAN(Model):
    def __init__(self, placeholders, input_dim, size_gene, latent_factor_num, **kwargs):
        super(ZINB_GRAN, self).__init__(**kwargs)
        self.inputs = placeholders['features']
        self.placeholders = placeholders
        self.size_gene = size_gene
        self.input_dim = input_dim
        self.latent_factor_num = latent_factor_num
        self.optimizer = tf.train.AdamOptimizer(learning_rate=FLAGS.learning_rate)
        self.layers = []
        self.build()

    def _loss(self):
        for var in self.layers[0].vars.values():
            self.loss += FLAGS.weight_decay * tf.nn.l2_loss(var)
        # self.loss += euclidean_loss(self.outputs, self.placeholders['labels'])
        self.loss += masked_accuracy(self.outputs, self.placeholders['labels'], 
                                    self.placeholders['labels_mask'],
                                    self.placeholders['negative_mask'])

    def _accuracy(self):
        # self.accuracy = euclidean_loss(self.outputs, self.placeholders['labels'])
        self.accuracy = masked_accuracy(self.outputs, self.placeholders['labels'], 
                                       self.placeholders['labels_mask'],
                                       self.placeholders['negative_mask'])

    def _build(self):
        # Encoder 1
        self.layers.append(Encoder(input_dim=self.input_dim,
                                   output_dim=FLAGS.hidden1,
                                   gene_size=self.size_gene,
                                   dropout=0.,
                                   featureless=False,
                                   placeholders=self.placeholders))

        # Encoder 2 (produces latent factors)
        self.layers.append(Encoder(input_dim=FLAGS.hidden1,
                                   output_dim=self.latent_factor_num,
                                   gene_size=self.size_gene,
                                   dropout=FLAGS.dropout,
                                   featureless=False,
                                   placeholders=self.placeholders,
                                   act=lambda x: x))
        
        # Decoder
        self.layers.append(Decoder(size1=self.size_gene,
                                   latent_factor_num=self.latent_factor_num,
                                   placeholders=self.placeholders))
        
        
        
    def predict(self):
        return self.outputs



class Discriminator(tf.keras.Model):
    def __init__(self, input_dim, hidden_dim1, hidden_dim2, output_dim):
        super(Discriminator, self).__init__()
        self.hidden_layer1 = tf.keras.layers.Dense(hidden_dim1, activation='relu')
        self.hidden_layer2 = tf.keras.layers.Dense(hidden_dim2, activation='relu')
        self.output_layer = tf.keras.layers.Dense(output_dim, activation='sigmoid')

    def build(self, input_shape):
        # Explicitly set input shape to trigger weight initialization
        self.hidden_layer1.build(input_shape)
        self.hidden_layer2.build((input_shape[0], self.hidden_layer1.units))
        self.output_layer.build((input_shape[0], self.hidden_layer2.units))

    def call(self, inputs):
        hidden1 = self.hidden_layer1(inputs)
        hidden2 = self.hidden_layer2(hidden1)
        output = self.output_layer(hidden2)
        return output





