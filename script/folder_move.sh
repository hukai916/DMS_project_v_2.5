# Clean the folders for pipeline analysis
dir=$PWD
workdir="$(dirname "$dir")"
outdir="$workdir"/result_example/"$1"
cd "$workdir"

mkdir "$outdir"
mv TemFolder "$outdir"
mv data_sra "$outdir"
mv enrich2_input "$outdir"
mv enrich2_json "$outdir"
mv first_mapper "$outdir"
mv merge "$outdir"
mv plot_input "$outdir"
mv plot_output "$outdir"
mv qc_INDEL "$outdir"
mv qc_library "$outdir"
mv second_mapper "$outdir"
mv ref "$outdir"

echo "Done moving!"
