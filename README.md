# Salience_Feature_Interactions
The main.m script demonstrates how the model learns optimized interaction weights across a set of acoustic features—envelope, pitch, low-frequency spectrogram, high-frequency spectrogram, scale, and rate. The framework first computes a set of surprisal values for each input feature, then applies logistic regression to learn feature interaction weights that best predict the presence of a deviant token within the audio trials.

The demonstration uses half of the trials from either the Music or Nature experiment for training, and the remaining half for evaluation. Finally, the script visualizes the learned weights and plots the ROC curve to assess the model's performance.

Please make sure to add https://github.com/JHU-LCAP/DREX-model and the rest of the folders in this repository to MATLAB's path before running this code. 
