pyinstaller --onefile --noconfirm -n apscale_blast-x64-macosx __main__.py
rm -rf build
tar -czf dist/apscale_blast-x64-macosx.zip dist/apscale_blast-x64-macosx
rm -rf dist/apscale_blast-x64-macosx
rm apscale_blast-x64-macosx.spec

# ./apscale_blast-x64-macosx blastn -database /Users/tillmacher/Desktop/_dev/PyInstaller_tests/APSCALE_blast/test_data/MIDORI2_UNIQ_NUC_GB259_srRNA_BLAST -query_fasta /Users/tillmacher/Desktop/_dev/PyInstaller_tests/APSCALE_blast/test_data/physalia_bf3_apscale_apscale_OTU_table_filtered.fasta
# ./apscale_blast-x64-macosx filter -database /Users/tillmacher/Desktop/_dev/PyInstaller_tests/APSCALE_blast/test_data/MIDORI2_UNIQ_NUC_GB259_srRNA_BLAST -blastn_folder /Users/tillmacher/Desktop/_dev/PyInstaller_tests/APSCALE_blast/dist/dist/blastn_physalia_bf3_apscale_apscale_OTU_table_filtered_07_01_24