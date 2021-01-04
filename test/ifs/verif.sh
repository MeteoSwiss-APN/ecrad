files=$(ls ecrad_meridian*_out.nc)
echo $files
for file in $files; do
  echo "Verifying ... " $file
  python3 -m recursive_diff.ncdiff $file refs/$file --rtol 1e-7
done

