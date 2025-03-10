function parse_yaml {    local prefix=$2;    local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034');    sed -ne "s|^\($s\):|\1|"         -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p"         -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |    awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'; }


for y in sample_*_rev/outs/*CB.yaml; do
  f=`dirname $y`
  echo $f
  eval $(parse_yaml $y)
  if [ -f "$f/cellbender_filtered.h5" ] ; then
    echo "Filtered results found; skipping"
  else
    cd $f
    cellbender remove-background --input raw_feature_bc_matrix.h5 --output cellbender_filtered.h5ad --cuda --expected-cells $expected_cells \
      --total-droplets-included $total_droplets_included --exclude-feature-types Peaks --low-count-threshold $low_count_threshold
    cd ../../
  fi
done
