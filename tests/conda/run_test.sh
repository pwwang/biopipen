if [[ -z "$1" ]]; then
    echo "Usage: $0 <proc_test_directory> [verbose (default: 0)]"
    echo "Example: $0 tests/test_misc/File2Proc 1"
    exit 1
fi

PROC_TEST_DIR=$1
PROCESS=$(basename $PROC_TEST_DIR)
VERBOSE=$2
if [[ -z "$VERBOSE" ]]; then
    VERBOSE=0
fi

echo "- Testing process: $PROCESS ..."

RUNFILE=$PROC_TEST_DIR/run.env
ENVNAME="base"
ARGS=""
LOCAL_ONLY="false"
if [[ -f $RUNFILE ]]; then
    source $RUNFILE
fi

if [[ "$LOCAL_ONLY" == "true" ]]; then
    echo "  Skipped"
    exit 0
else
    CMD="conda run --no-capture-output -n $ENVNAME poetry run python $PROC_TEST_DIR/test.py $ARGS"
    echo "  Running: $CMD"

    if [[ $VERBOSE -eq 1 ]]; then
        conda run --no-capture-output -n $ENVNAME poetry run python $PROC_TEST_DIR/test.py $ARGS
    else
        conda run --no-capture-output -n $ENVNAME poetry run python $PROC_TEST_DIR/test.py $ARGS > /dev/null 2>&1
    fi

    if [[ $? -eq 0 ]]; then
        echo "  v Success :)"
    else
        echo "  x Failure :("
        exit 1
    fi
fi
