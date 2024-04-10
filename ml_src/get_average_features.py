import ml_src.FEAT as FEAT

def get_average_features(sequence:str)->list:
	"""
	Method: Fills feature columns in a dataframe with net averages

	Input:
	  - sequence: amino acid string
	
	Output:
    - list of averages
	"""

	# get net cumulative sums
	net_cumulative_sum = min(len(sequence), 900)

	gravy = hydrophobicity = exposed = disorder = bulkiness = interface = 0.0
	for i, aa in enumerate(sequence):
		if i == net_cumulative_sum: break

		if aa.upper() in FEAT.INTERFACE_DIC:
			gravy += FEAT.GRAVY_DIC[aa.upper()]
			hydrophobicity += FEAT.HYDRO_DIC[aa.upper()]
			exposed += FEAT.EXPOSED_DIC[aa.upper()]
			disorder += FEAT.DISORDER_DIC[aa.upper()]
			bulkiness += FEAT.BULKY_DIC[aa.upper()]
			interface += FEAT.INTERFACE_DIC[aa.upper()]

	# return averages
	return [gravy           / net_cumulative_sum,
          hydrophobicity  / net_cumulative_sum,
          exposed 		    / net_cumulative_sum,
          disorder 		    / net_cumulative_sum,
          bulkiness 		  / net_cumulative_sum,
          interface 		  / net_cumulative_sum]