Data code for the manuscript 'Rapid threat detection: temporal priority in processing nociceptive intensity over location across species'.
Code_for_for_EEG_laser_dataset contains the tensor decomposition code in Experiment 4.
./Code_for_IntensityPriority.m: the code for constructing the spectral data from raw EEG data before following tensor composition.
./TCA_eval.py and ./TCA_shuffle.py: the code for tensor decomposition on spectral data and shuffled data using tensortool in python
./Tensor_decomposition_for_eeg.m: code for tensor classification, statistical comparison and figure plotting

Code_for_for_EEG_two_modality contains the tensor decomposition code in Experiment 5.
./Code_for_IntensityPoriority.m ./TCA_eval.py and ./TCA_shuffle.py : similar codes in Code_for_IntensityPriority
./Tensor_decomposition_for_EEG_twomodality.m: code for tensor classification, statistical comparison and figure plotting

Code_for_spike contains the tensor decomposition code and neural trajectory analysis in Experiment 6
./TCA_cal_silicon_twomodality.py and ./TCA_shuffle_silicon_twomodality the code for tensor decomposition on structed SPKdata (.mat) and shuffled data using tensortool in python.
./TCAdecomposition_for_silicon_twomodality.m  code for tensor classification, statistical comparison and figure plotting on SPKdata.
./Neuraltrajectory.m neural trajectory analysis to show the encoding preference of SPKdata on nociceptive intensity and location.
