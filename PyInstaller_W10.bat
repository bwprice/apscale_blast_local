pyinstaller --onefile --noconfirm -n apscale_blast-x64-windows __main__.py
rmdir /Q /S build
tar -a -c -f dist\\apscale_blast-x64-windows.zip dist\\apscale_blast-x64-windows
rmdir /Q /S dist\\apscale_blast-x64-windows
del /f apscale_blast-x64-windows.spec
