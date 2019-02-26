#![feature(type_alias_enum_variants)]


mod csidh;

mod secret_sharing;


fn main() {
    secret_sharing::secret_share_demo(1);
}