# Import necessary libraries
import numpy as np
import tensorflow.compat.v1 as tf
import pandas as pd
tf.disable_v2_behavior()

index = str(sys.argv[1])

# Get the directory of the current Python file
current_directory = os.path.dirname(os.path.abspath(__file__))


# Define the true confusion matrix
true_confusion_matrix = np.array([
    [0.9, 0.05, 0.05],
    [0.05, 0.9, 0.05],
    [0.05, 0.05, 0.9]
])


# Step 1: Load real data
covariates_file = os.path.join(current_directory, f'trainX{index}.csv')
noisy_labels_file = os.path.join(current_directory, f'trainY{index}.csv')

# Load data
covariates_train = pd.read_csv(covariates_file).values.astype(np.float32)
noisy_labels_train = pd.read_csv(noisy_labels_file).values.astype(np.int32)

# Ensure data dimensions
num_samples_train, vector_size = covariates_train.shape
num_annotators = noisy_labels_train.shape[1]

# Convert noisy labels to one-hot encoding
num_classes = len(np.unique(noisy_labels_train))
annotator_labels_train = np.zeros((num_samples_train, num_annotators, num_classes), dtype=np.float32)
for i in range(num_samples_train):
    for j in range(num_annotators):
        annotator_labels_train[i, j, :] = np.eye(num_classes)[noisy_labels_train[i, j]]



# Step 1: Load real data
covariates_file = os.path.join(current_directory, f'testX{index}.csv')
noisy_labels_file = os.path.join(current_directory, f'testY{index}.csv')

# Load data
covariates_test = pd.read_csv(covariates_file).values.astype(np.float32)
noisy_labels_test = pd.read_csv(noisy_labels_file).values.astype(np.int32)

# Ensure data dimensions
num_samples_test= covariates_test.shape[0]

# Convert noisy labels to one-hot encoding
annotator_labels_test = np.zeros((num_samples_test, num_annotators, num_classes), dtype=np.float32)
for i in range(num_samples_test):
    for j in range(num_annotators):
        annotator_labels_test[i, j, :] = np.eye(num_classes)[noisy_labels_test[i, j]]


# Define scale values
scale_values = np.arange(0.20, 0.80, 0.05)


# Placeholders
vector_input = tf.placeholder(tf.float32, shape=(None, vector_size), name="vector_input")
annotator_labels_input = tf.placeholder(
    tf.float32, shape=(None, num_annotators, num_classes), name="annotator_labels_input"
)

# Classifier
def classifier(vectors):
    hidden1 = tf.keras.layers.Dense(64, activation=tf.nn.relu, name="hidden1")(vectors)
    hidden2 = tf.keras.layers.Dense(32, activation=tf.nn.relu, name="hidden2")(hidden1)
    logits = tf.keras.layers.Dense(num_classes, name="output")(hidden2)
    return logits

logits = classifier(vector_input)

# Confusion matrix estimators
def confusion_matrix_estimators(num_annotators, num_classes):
    w_init = tf.constant(
        np.stack([6.0 * np.eye(num_classes) - 5.0 for _ in range(num_annotators)]),
        dtype=tf.float32,
    )
    rho = tf.Variable(w_init, name="rho")
    rho = tf.nn.softplus(rho)
    confusion_matrices = tf.divide(rho, tf.reduce_sum(rho, axis=-1, keepdims=True))
    return confusion_matrices

confusion_matrices = confusion_matrix_estimators(num_annotators, num_classes)

# Loss function
def cross_entropy_over_annotators(labels, logits, confusion_matrices):
    labels = tf.cast(labels, dtype=tf.float32)
    losses_all_annotators = []
    for idx, labels_annotator in enumerate(tf.unstack(labels, axis=1)):
        preds_true = tf.nn.softmax(logits)
        preds_annotator = tf.matmul(preds_true, confusion_matrices[idx, :, :])
        preds_clipped = tf.clip_by_value(preds_annotator, 1e-10, 0.9999999)
        loss = tf.reduce_sum(-labels_annotator * tf.log(preds_clipped), axis=-1)
        losses_all_annotators.append(loss)
    losses_all_annotators = tf.stack(losses_all_annotators, axis=1)
    has_labels = tf.reduce_sum(labels, axis=2)
    losses_all_annotators = losses_all_annotators * has_labels
    return tf.reduce_mean(tf.reduce_sum(losses_all_annotators, axis=1))

weighted_cross_entropy = cross_entropy_over_annotators(annotator_labels_input, logits, confusion_matrices)
trace_norm = tf.reduce_mean(tf.trace(confusion_matrices))

# Optimizer
scale = tf.placeholder(tf.float32, name="scale")
total_loss = weighted_cross_entropy + scale * trace_norm
optimizer = tf.train.AdamOptimizer(learning_rate=0.001)
train_op = optimizer.minimize(total_loss)


# Training loop for each scale value
results = []
num_epochs = 1000
batch_size = 4
true_labels_train = np.array([0] * 100 + [1] * 100 + [2] * 100)
true_labels_test = np.array([0] * 30 + [1] * 30 + [2] * 30)

with tf.Session() as sess:
    for scale_value in scale_values:
        sess.run(tf.global_variables_initializer())
        

        for epoch in range(num_epochs):
            for i in range(0, num_samples_train, batch_size):
                batch_vectors = covariates_train[i:i + batch_size]
                batch_annotator_labels = annotator_labels_train[i:i + batch_size]
                sess.run(train_op, feed_dict={
                    vector_input: batch_vectors,
                    annotator_labels_input: batch_annotator_labels,
                    scale: scale_value
                })

        # Evaluate accuracy for training
        logits_output_train = sess.run(logits, feed_dict={vector_input: covariates_train})
        predicted_classes_train = np.argmax(logits_output_train, axis=1)
        accuracy_train = np.mean(predicted_classes_train == true_labels_train)

        # Evaluate accuracy for testing
        logits_output_test = sess.run(logits, feed_dict={vector_input: covariates_test})
        predicted_classes_test = np.argmax(logits_output_test, axis=1)
        accuracy_test = np.mean(predicted_classes_test == true_labels_test)


        # Calculate metrics
        estimated_confusion_matrices = sess.run(confusion_matrices)
        f_norms = [np.linalg.norm(cm - true_confusion_matrix, ord='fro') for cm in estimated_confusion_matrices]
        average_f_norm = np.mean(f_norms)
        traces = np.trace(estimated_confusion_matrices, axis1=1, axis2=2)
        average_trace = np.mean(traces)
        final_loss = sess.run(total_loss, feed_dict={
            vector_input: covariates_train,
            annotator_labels_input: annotator_labels_train,
            scale: scale_value
        })

        # Store results
        results.append((scale_value, average_f_norm, average_trace, final_loss, accuracy_train, accuracy_test))



# Save `results`
results_df = pd.DataFrame(results, columns=['scale_value', 'average_f_norm', 'average_trace', 'final_loss', 'accuracy_train', 'accuracy_test'])
results_path = os.path.join(current_directory, f'results{index}.csv')
results_df.to_csv(results_path, index=False)