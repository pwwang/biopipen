# Keras Model Migration Guide

## Overview

The original `TrainedEncoder.h5` model was created with Keras <v3. This directory now contains updated code that supports both Keras v2 and Keras v3.

## What Changed

### 1. Updated `BriseisEncoder.py`

The script now:
- Automatically detects the Keras version
- Supports loading both `.h5` (legacy) and `.keras` (Keras v3) format models
- First checks for a `.keras` version of the model when running with Keras v3
- Falls back to legacy loading methods if needed

### 2. New Migration Script

A new script `migrate_keras_model.py` is provided to convert the old `.h5` model to the new `.keras` format.

## Migration Instructions

### Option 1: Automatic (Recommended)

The `BriseisEncoder.py` script will automatically handle loading the model:

1. If using Keras v3 and a `.keras` version exists, it will use that
2. If using Keras v3 with only a `.h5` file, it will attempt legacy loading
3. If using Keras v2, it will load normally

**No action required** - the script handles everything automatically.

### Option 2: Manual Migration (Optional)

If you want to pre-convert the model to Keras v3 format:

1. Ensure you have Keras v3 installed:
   ```bash
   pip install --upgrade keras>=3.0
   ```

2. Run the migration script:
   ```bash
   python migrate_keras_model.py -i TrainedEncoder.h5 -o TrainedEncoder.keras
   ```

3. The script will create `TrainedEncoder.keras` which will be automatically detected and used

### Option 3: Using the R Script

When using the TESSA R script (`TESSA.R`), you can optionally specify the model path:

The R script will pass the model path to `BriseisEncoder.py`, which will handle the model format automatically.

## Compatibility

- **Keras v2**: Works with `.h5` files (original format)
- **Keras v3**: Works with both `.h5` (legacy) and `.keras` (native) files
  - `.keras` files load faster and are the preferred format for Keras v3
  - `.h5` files are automatically converted on-the-fly using legacy loaders

## Troubleshooting

### Error: "Failed to load model"

If you encounter errors loading the model:

1. Check your Keras version:
   ```python
   import keras
   print(keras.__version__)
   ```

2. Try the manual migration (Option 2 above)

3. If issues persist, ensure you have compatible dependencies:
   ```bash
   pip install --upgrade tensorflow keras numpy pandas
   ```

### Error: "Argument(s) not recognized: {'lr': ...}"

This error occurs when Keras v3 tries to load an optimizer saved with Keras v2 parameter names (e.g., `lr` instead of `learning_rate`). The BriseisEncoder.py script handles this automatically by loading models without compilation (`compile=False`), which avoids optimizer parameter conflicts. The model will still work correctly for inference.

If you need the compiled model for training, you'll need to:
1. Migrate the model to `.keras` format using the migration script
2. Recompile the model manually with the correct Keras v3 optimizer parameters

### Performance Notes

- Keras v3 `.keras` format models load faster than legacy `.h5` files
- The first-time legacy load may be slower but subsequent runs will use cached data
- For production use, we recommend pre-converting to `.keras` format using the migration script

## Files

- `TrainedEncoder.h5` - Original Keras v2 model (legacy format)
- `TrainedEncoder.keras` - Keras v3 model (created after migration)
- `migrate_keras_model.py` - Script to convert between formats
- `BriseisEncoder.py` - Updated encoder script with automatic version detection
- `KERAS_MIGRATION.md` - This file

## Additional Resources

- [Keras 3 Migration Guide](https://keras.io/guides/migrating_to_keras_3/)
- [Keras Saving and Serialization](https://keras.io/guides/serialization_and_saving/)
