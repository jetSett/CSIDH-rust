#![feature(type_alias_enum_variants)]
#![feature(specialization)]


mod csidh;

mod secret_sharing;
mod ot;

fn main() {
    // ot::oblivious_transfert_demo(30, [42, 0]);
    ot::oblivious_transfert_fake(30, [42, 0]);
    //secret_sharing::secret_share_demo(100);
}