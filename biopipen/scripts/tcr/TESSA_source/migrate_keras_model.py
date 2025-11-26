#!/usr/bin/env python3
"""
Migrate Keras <v3 model (TrainedEncoder.h5) to Keras v3 format.

This script converts the old Keras model to the new .keras format compatible with Keras v3.
"""

import os
import sys
import argparse


def migrate_model(old_model_path, new_model_path):
    """
    Migrate a Keras <v3 .h5 model to Keras v3 .keras format.

    Args:
        old_model_path: Path to the old .h5 model file
        new_model_path: Path to save the new .keras model file
    """
    try:
        # Try importing keras v3
        import keras

        # Check if we're using Keras 3
        keras_version = keras.__version__
        print(f"Using Keras version: {keras_version}")

        if keras_version.startswith('3.'):
            print("Keras v3 detected. Attempting to load legacy model...")

            # In Keras 3, we need to use the legacy loader for old h5 files
            try:
                from keras.src.legacy.saving import legacy_h5_format
                model = legacy_h5_format.load_model_from_hdf5(old_model_path)
                print(f"Successfully loaded legacy model from {old_model_path}")
            except (ImportError, AttributeError):
                # Fallback: try using keras.saving.load_model with compile=False
                print("Attempting to load with keras.saving.load_model...")
                model = keras.saving.load_model(old_model_path, compile=False)
                print(f"Successfully loaded model from {old_model_path}")

            # Save in the new Keras 3 format
            model.save(new_model_path)
            print(f"Successfully saved model to {new_model_path} in Keras v3 format")

            # Verify the new model can be loaded
            print("Verifying new model...")
            test_model = keras.saving.load_model(new_model_path)
            print("Verification successful!")

            return True
        else:
            print(f"ERROR: This script requires Keras v3, but found version {keras_version}")
            print("Please upgrade keras: pip install --upgrade keras")
            return False

    except ImportError as e:
        print(f"ERROR: Could not import keras: {e}")
        print("Please install keras: pip install keras")
        return False
    except Exception as e:
        print(f"ERROR: Failed to migrate model: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Migrate Keras <v3 model to Keras v3 format"
    )
    parser.add_argument(
        "-i", "--input",
        default="TrainedEncoder.h5",
        help="Input .h5 model file (default: TrainedEncoder.h5)"
    )
    parser.add_argument(
        "-o", "--output",
        default="TrainedEncoder.keras",
        help="Output .keras model file (default: TrainedEncoder.keras)"
    )

    args = parser.parse_args()

    # Get the directory of this script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Resolve paths
    if not os.path.isabs(args.input):
        input_path = os.path.join(script_dir, args.input)
    else:
        input_path = args.input

    if not os.path.isabs(args.output):
        output_path = os.path.join(script_dir, args.output)
    else:
        output_path = args.output

    # Check if input file exists
    if not os.path.exists(input_path):
        print(f"ERROR: Input file not found: {input_path}")
        sys.exit(1)

    print(f"Input model: {input_path}")
    print(f"Output model: {output_path}")

    # Migrate the model
    success = migrate_model(input_path, output_path)

    if success:
        print("\n" + "="*60)
        print("Migration completed successfully!")
        print("="*60)
        sys.exit(0)
    else:
        print("\n" + "="*60)
        print("Migration failed!")
        print("="*60)
        sys.exit(1)


if __name__ == "__main__":
    main()
