import sys
import os
from argparse import ArgumentParser
import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(SCRIPT_DIR)

from ml_src.ml_model import ML_Model
from gen_src.fasta_content import Fasta_Content
from ml_src.get_average_features import get_average_features

## take in:
##    1) secreted proteins fasta file
##    2) (optional) model file path (will default to best Random Forest)

## RUN LIKE THIS:
##    Using conda: python3 effectorO.py -i {INPUT_FASTA_PATH}

## output:
##    1) csv conatining (by column):
##			 - protein_id from FASTA file
##			 - sequence from FASTA file
##			 - prediction - 0=non_effector | 1=effector
##			 - probability that a sequence is an effector
##			 - meaning - verbose string version of prediction
##    2) fasta file of predicted effectors


def __get_seq_features(fasta_content:Fasta_Content):
	print("Calculating sequence features...")
	seq_features = np.array([get_average_features(seq) for seq in fasta_content.get_sequences()])

	# Ensure seq_features is a 2D array
	if seq_features.ndim == 1:
			seq_features = np.array(seq_features).reshape(1, -1)
	elif seq_features.ndim == 0:
		raise ValueError("Getting average features did not produce a result")
	elif seq_features.ndim > 2:
			raise ValueError("Sequence features should be represented in either a 1D or 2D array")

	return seq_features

def _analyze_fasta_cli(fasta_input:str):
	print("Parsing FASTA file...")
	fasta_content = Fasta_Content()
	fasta_content.parse_fasta_file(fasta_input)
	seq_features = __get_seq_features(fasta_content)
	return (fasta_content, seq_features)

def _analyze_fasta_api(filename:str, content:str):
	fasta_content = Fasta_Content()
	fasta_content.parse_fasta_content(filename, content)
	seq_features = __get_seq_features(fasta_content)
	return (fasta_content, seq_features)

def _create_output_dir(output_dir:str=None):
	print("Creating a directory of EffectorO outputs...")
	if not output_dir: output_dir = "effectoro_results"
	if os.path.exists(output_dir):
		print("Warning: directory 'effectoro_results' already exists. Placing contents into the existing directory...")
	os.makedirs(output_dir, exist_ok=True)
	return output_dir

def _create_csv(predictions_df:pd.DataFrame, filepath:str):
	print("Constructing a CSV containing information on sequences from the FASTA file...")
	predictions_df.to_csv(filepath, index=False)

def _create_fasta(predictions_df:pd.DataFrame, filepath:str, effector_prediction_type:str=None):
	print("Constructing a FASTA file with just predicted effectors...")
	with open(filepath, 'w') as outfile:
		predicted_effectors = (
			(predictions_df[predictions_df['meaning'].astype(str) == effector_prediction_type])
				if effector_prediction_type else predictions_df
		)
		seqs_to_write = [
				f">{row['protein_id']} {row['meaning']} probability={row['probability']}\n{row['sequence']}\n"
				for _, row in predicted_effectors.iterrows()
		]
		outfile.writelines(seqs_to_write)

def _run_effectoro_ml(fasta_content_to_predict, seq_features, model_path=None):
	model = ML_Model(model_path, SCRIPT_DIR)

	print("Sequences to run secreted oomycete-trained " +
				f"{'Random Forest' if (model.model_name == 'RF_88_best') else (model.model_name)} " +
				f"effector classifier on: {len(seq_features)}")
	if model.model_name == 'RF_88_best':
		print("\n**NOTES**:\n\n" +
					"Positive training dataset: " +
					"~100 experimentally validated oomycete avirulence effectors\n" +
					"Negative training dataset: " +
					"~100 secreted orthologous oomycete genes\n" +
					"\n**END OF NOTES**\n")

	predictions_df = model.predict_effectors(seq_features, fasta_content_to_predict)
	return predictions_df

def run_effectoro_ml_through_cli(input_fasta:str, model_path, output_dir:str):
	fasta_content_to_predict, seq_features = _analyze_fasta_cli(input_fasta)
	predictions = _run_effectoro_ml(fasta_content_to_predict, seq_features, model_path)

	output_dir = _create_output_dir(output_dir)
	csv_filepath = os.path.join(output_dir, fasta_content_to_predict.get_fasta_filename() + ".effector_classification_table.csv")
	_create_csv(predictions, csv_filepath)
	pef_fasta_filepath = os.path.join(output_dir, fasta_content_to_predict.get_fasta_filename() + ".predicted_effectors.fasta")
	_create_fasta(predictions, pef_fasta_filepath, "predicted_effector")
	pnef_fasta_filepath = os.path.join(output_dir, fasta_content_to_predict.get_fasta_filename() + ".predicted_non-effectors.fasta")
	_create_fasta(predictions, pnef_fasta_filepath, "predicted_non-effector")
	all_fasta_filepath = os.path.join(output_dir, fasta_content_to_predict.get_fasta_filename() + ".all_effector_predictions.fasta")
	_create_fasta(predictions, all_fasta_filepath)

	print(f"\nOutput fasta of predicted effectors available in:\n{pef_fasta_filepath}\n" +
			  f"\nOutput fasta of predicted non-effectors available in:\n{pnef_fasta_filepath}\n" +
				f"\nOutput fasta of all effector predictions available in:\n{all_fasta_filepath}\n" +
				f"\nDetailed CSV with seqeunce IDs | sequences | predictions | probabilities in:\n{csv_filepath}\n")

def run_effectoro_ml_through_api(file_name: str, fasta_content:str)->pd.DataFrame:
	fasta_content_to_predict, seq_features = _analyze_fasta_api(file_name, fasta_content)
	predictions = _run_effectoro_ml(fasta_content_to_predict, seq_features)
	return predictions


if __name__ == "__main__":
	parser = ArgumentParser(prog="A Python script to run EffectorO.")
	parser.add_argument("-i", "--input_fasta", type=str, help="Input FASTA file", required=True)
	parser.add_argument("-m", "--model_path", type=str, help="Path to trained model to be used [default: provided Random Forest model]", required=False)
	parser.add_argument("-o", "--output_dir", type=str, help="Output directory")
	parser.add_argument("-s", "--suppress_prints", action="store_true", help="Suppress all print statements")
	args = parser.parse_args()

	# suppress print statements if indicated
	if args.suppress_prints:
		sys.stdout = open(os.devnull, 'w')
		sys.stderr = open(os.devnull, 'w')

	# run the effectoro-ML script
	run_effectoro_ml_through_cli(args.input_fasta, args.model_path, args.output_dir)
