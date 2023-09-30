if [[ -z "$1" ]]; then
    echo "Usage: $0 <proc_test_directory> [verbose (default: 0)]"
    echo "Example: $0 tests/test_misc/File2Proc 1"
    exit 1
fi

PROC_TEST_DIR=$1
PROCESS=$(basename "$PROC_TEST_DIR")
VERBOSE_OR_FORCE=$2
FORCE_OR_VERBOSE=$3
VERBOSE="false"
FORCE="false"
if [[ "$VERBOSE_OR_FORCE" == "VERBOSE="* ]]; then
    # remove the "VERBOSE=" prefix
    VERBOSE=${VERBOSE_OR_FORCE:8}
fi
if [[ "$VERBOSE_OR_FORCE" == "FORCE="* ]]; then
    # remove the "FORCE=" prefix
    FORCE=${VERBOSE_OR_FORCE:6}
fi
if [[ "$FORCE_OR_VERBOSE" == "VERBOSE="* ]]; then
    # remove the "VERBOSE=" prefix
    VERBOSE=${FORCE_OR_VERBOSE:8}
fi
if [[ "$FORCE_OR_VERBOSE" == "FORCE="* ]]; then
    # remove the "FORCE=" prefix
    FORCE=${FORCE_OR_VERBOSE:6}
fi
if [[ -z "$VERBOSE" ]]; then
    VERBOSE="false"
fi
if [[ -z "$FORCE" ]]; then
    FORCE="false"
fi

echo "- Testing process: $PROCESS ..."

RUNFILE="$PROC_TEST_DIR/run.env"
ENVNAME="base"
ARGS=""
LOCAL_ONLY="false"
if [[ -f $RUNFILE ]]; then
    # shellcheck disable=SC1090
    source "$RUNFILE"
fi

if [[ "$LOCAL_ONLY-$FORCE" == "true-false" ]]; then
    echo "  Skipped"
    exit 0
else
    CMD="conda run --no-capture-output -n $ENVNAME poetry run python $PROC_TEST_DIR/test.py $ARGS"
    echo "  Running: $CMD"

    if [[ $VERBOSE == "true" ]]; then
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
