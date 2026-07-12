use fs_gum_topo::hopf::hopf_whitehead;
use fs_gum_topo::referee::hopfion_director;
use fs_gum_topo::Bc;

fn main() {
    for n in [32usize, 64] {
        let f = hopfion_director(n, 2.4, Bc::Periodic);
        match hopf_whitehead(&f) {
            Ok(w) => println!("N={n}: H={:.9} raw={:.9} mean_b={:?} max_b={:.3e} curl={:.3e} im={:.3e}", w.hopf, w.raw, w.mean_b, w.max_b, w.curl_rel_err, w.max_imag),
            Err(e) => println!("N={n}: ERR {e}"),
        }
    }
}
