# Clean the folders for pipeline analysis
dir=$PWD
workdir="$(dirname "$dir")"
cd "$workdir"

rm -r TemFolder data_sra enrich2_input enrich2_json enrich2_output first_mapper merge plot_input plot_output qc_INDEL qc_library second_mapper ref

echo "Done cleaning!"
