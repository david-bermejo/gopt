Write-Host "Building..."
Set-Location build
& "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Tools\Llvm\bin\clang++.exe" --std=c++2a -O3 ../main.cpp
Set-Location ..
Write-Host "Done."
Write-Host ""
Write-Host "Executing program..."
Write-Host ""
./build/a.exe