# Macros for ML training in D2H
## Requirements
In order to execute the training macro in this folder, the following python libraries are needed:
- [hipe4ml](https://github.com/hipe4ml/hipe4ml)
- [hipe4ml_converter](https://github.com/hipe4ml/hipe4ml_converter)

All configurables needed for sample preparation and ML training are respectively embedded in configuration files such as:
- `config_preparation_DplusToPiKPi.yml`;
- `config_training_DplusToPiKPi.yml`.

## Training
### Download training samples from Hyperloop
Download the derived `AO2D.root` (produced via a treeCreator task for example) and save it in a dedicated directory (one directory per workflow). In each directory, the output of the data preparation will be stored. A posteriori, inside the training configuration file, one shall specify the path to these directories. There is then no need to specify the files names as these ones are generated automatically (according to the analysed decay channel).

### Prepare samples

In order to prepare the samples the following script can be used while being in O2(Physics) environment:
```
python3 prepare_samples.py config_preparation_DplusToPiKPi.yml
```

Preselections can be applied at this level. This macro shall produce `.root` files containing labeled data samples for each input `AO2D` file.


### Perform training
In order to perform the training and produce the BDT models to be used in D2H, the following script can be used:
```
python3 train_models.py config_training_DplusToPiKPi.yml
```
Given the output directory set in `config_training_DplusToPiKPi.yml`, a directory is created for each pT bin (i.e. each model trained) and filled with:
- plots at data preparation level: variables distributions and correlations
- plots at training-testing level: BDT output scores, ROC curves, precision recall, features importance
- trained models in files containing `ModelHandler` prefix (in `pickle` and `onnx` formats). *The `.onnx` can be used for ML inference in O2Physics selection task.*
- model applied to test set in file containing `ModelApplied` suffix
