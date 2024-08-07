import os
import numpy as np
import pandas as pd
import joblib

from gen_src.fasta_content import Fasta_Content

class ML_Model:

  def __init__(self, model_path:str=None, cwd:str=""):
    print("Importing ML model...")
    model_path = str(model_path) if (model_path) \
      else os.path.join(cwd, "trained_models", "RF_88_best.sav")

    self.model_name = os.path.splitext(os.path.basename(model_path))[0]
    self.model = None
    with open(model_path, 'rb') as fmodel:
      self.model = joblib.load(fmodel)

  def predict_effectors(self, seq_features:np.ndarray, fasta_content:Fasta_Content):
    print("Predicting effectors amongst FASTA sequences...")

    # perform predictions
    predictions = self.model.predict(seq_features)
    probabilities = self.model.predict_proba(seq_features)

    # map predictions to meanings
    PREDICTION_MAP = {'0': "predicted_non-effector", '1': "predicted_effector"}
    meanings = np.array([PREDICTION_MAP[str(pred)] for pred in predictions])

    # create dataframe
    result_df = pd.DataFrame({
        "protein_id": fasta_content.get_ids(),
        "sequence": fasta_content.get_sequences(),
        "prediction": predictions,
        "probability": probabilities[:, 1],
        "meaning": meanings
    })

		# define dataframe column types
    result_df["protein_id"] = result_df["protein_id"].astype(str)
    result_df["sequence"] = result_df["sequence"].astype(str)
    result_df["prediction"] = result_df["prediction"].astype(int)
    result_df["probability"] = result_df["probability"].astype(float)
    result_df["meaning"] = result_df["meaning"].astype(str)

    # round probabilities
    result_df['probability'] = round(result_df['probability'], 2)

    # print counts of predicted classes
    print(f"\nCounts of predicted classes:\n{result_df['meaning'].value_counts().to_string(header=False)}\n")

    return result_df
