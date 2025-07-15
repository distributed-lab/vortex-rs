use std::path::PathBuf;

#[cfg(all(
    feature = "nightly-features",
    target_arch = "x86_64",
    target_feature = "avx512"
))]
fn main() {
    // Tell rustc where the Go‑built library lives
    let lib_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("native")
        .join("libs");

    println!("cargo:rustc-link-search=native={}", lib_dir.display());

    // Link against libsis_amd64.so  (if you have libsis_amd64.a, change to static)
    println!("cargo:rustc-link-lib=dylib=sis_amd64");

    // Linux convenience: embed rpath so tests/binaries find the .so at runtime
    #[cfg(target_os = "linux")]
    println!("cargo:rustc-link-arg=-Wl,-rpath,$ORIGIN/../native/lib");
}

#[cfg(all(not(target_feature = "avx512")))]
fn main() {}
