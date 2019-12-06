CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate py3-dms

MODE=single
SCALE=max
CONFIG=NULL
HELP=NULL

for i in "$@"
do
case $i in
    -m=*|--mode=*)
    MODE="${i#*=}"
    ;;
    -s=*|--scale=*)
    SCALE="${i#*=}"
    ;;
    -c=*|--config=*)
    CONFIG="${i#*=}"
    ;;
    -h*|--help*)
    HELP=True
    ;;
esac
done

if [ $HELP = True ]
then
python wrapper.py
exit 1
fi

if [ $CONFIG = NULL ]
then
echo "Missing/wrong parameters, check usage:
bash wrapper.sh -h."
exit 1
fi

if [ $MODE = single ] && [ $SCALE = max ]
then
python wrapper.py -c $CONFIG -m single -s max
elif [ $MODE = double ] && [ $SCALE = max ]
then
python wrapper.py -c $CONFIG -m double -s max
elif [ $MODE = single ] && [ $SCALE != max ]
then
python wrapper.py -c $CONFIG -m single -s $SCALE
elif [ $MODE = double ] && [ $SCALE != max ]
then
python wrapper.py -c $CONFIG -m double -s $SCALE
else
echo "Missing/wrong parameters, check usage:
bash wrapper.sh -h."
fi
