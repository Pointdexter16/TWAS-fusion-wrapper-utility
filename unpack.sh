
source $(conda info --base)/etc/profile.d/conda.sh && conda activate && pip install gdown && gdown --id 1wX-6wIR92a4SxgBS58mEknKcFG0HuhFo -O master.zip
unzip master.zip

wget https://data.broadinstitute.org/alkesgroup/FUSION/SUM/PGC2.SCZ.sumstats

mkdir WEIGHTS
cd WEIGHTS
wget https://data.broadinstitute.org/alkesgroup/FUSION/WGT/GTEx.Whole_Blood.tar.bz2
tar xjf GTEx.Whole_Blood.tar.bz2
rm GTEx.Whole_Blood.tar.bz2
cd ..
conda env create -f environmentLdsc.yml
conda env create -f environmentFiz.yml

echo "alias fuse='make -f $(pwd)/Makefile base=$(pwd)'" >> ~/.zshrc
source $(conda info --base)/etc/profile.d/conda.sh && conda activate ldsc && conda install -c conda-forge r-rcpp r-rcppeigen -y && Rscript inPlink.R
