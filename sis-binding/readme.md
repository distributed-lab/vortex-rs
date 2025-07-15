# RingSIS bindings for Rust

To compile:
```shell
docker run --rm -v "$PWD":/src -w /src --platform linux/amd64 -it --entrypoint bash golang:1.24-bullseye
```

Then
```shell
export GOAMD64='v4'
go env -w GOAMD64='v4'
export CC=x86_64-linux-gnu-gcc
export CXX=x86_64-linux-gnu-g++
go build -buildmode=c-shared -o /src/sis_amd64.so ./binding
```