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

CONTAINER_FILE=$PROC_TEST_DIR/container.txt
if [[ ! -f $CONTAINER_FILE ]]; then
    echo "  No container.txt found, using host to test"
    cmd="python $PROC_TEST_DIR/pipeline.py"
    echo "  Running: $cmd"
else
    container="biopipen/$(cat $CONTAINER_FILE)"
    echo "  Using container: $container"
    cmd="docker run --rm -v $PWD:/workdir -w /workdir --entrypoint bash $container entrypoint.sh $PROC_TEST_DIR/pipeline.py"
    echo "  Running: $cmd"
fi

if [[ $VERBOSE -eq 1 ]]; then
    eval $cmd
else
    eval $cmd > /dev/null 2>&1
fi

if [[ $? -eq 0 ]]; then
    echo "  v Success :)"
else
    echo "  x Failed :("
    exit 1
fi
