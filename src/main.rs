#![feature(type_alias_enum_variants)]


mod csidh;

mod secret_sharing;
mod ot;
use elliptic_curve_algorithms::field::Field;
use ot::*;

fn main() {
    ot::oblivious_transfert_demo(10, [42, 0]);
    //secret_sharing::secret_share_demo(1);
}