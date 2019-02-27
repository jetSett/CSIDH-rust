use crate::csidh::*;

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use elliptic_curve_algorithms::finite_fields::*;

elliptic_curve_algorithms::declare_finite_field!(K, Integer, 
            Integer::from_str_radix("14841476269619", 10).unwrap(), 
            _m);

type PublicKey = K;
type SecretKey = Vec<i32>;

fn hash_k(val : K) -> u64{
    let mut hasher = DefaultHasher::new();
    Hash::hash_slice(&[val.clone()], &mut hasher);
    hasher.finish()
}


struct Alice<'a>{
    pub pk0 : PublicKey,
    sk0 : SecretKey,

    pub pk1 : PublicKey,
    sk1 : SecretKey,

    inst : &'a CSIDHInstance<K>
}

impl<'a> Alice<'a> {
    fn step_1(inst : &'a CSIDHInstance<K>, gen_param : i32) -> Alice<'a>{
        let k0 = inst.sample_keys(gen_param);
        let k1 = inst.sample_keys(gen_param);

        Alice{
            pk0: k0.0,
            sk0: k0.1,

            pk1: k1.0,
            sk1: k1.1,

            inst
        }
    }

    fn step_2(&self, s : [u64; 2], pk_bob : PublicKey) -> [u64; 2]{
        let mut sk0_inv = self.sk0.clone();
        let mut sk1_inv = self.sk1.clone();

        for i in 0..sk0_inv.len(){
            sk0_inv[i] = -sk0_inv[i];
        }

        for i in 0..sk1_inv.len(){
            sk1_inv[i] = -sk1_inv[i];
        }

        [s[0]^hash_k(self.inst.class_group_action(pk_bob.clone(), sk0_inv)),
        s[1]^hash_k(self.inst.class_group_action(pk_bob.clone(), sk1_inv))]
    }

}

struct Bob<'a> {
    pk: PublicKey,
    sk: SecretKey,

    k: usize,

    inst : &'a CSIDHInstance<K>
}


impl<'a> Bob<'a> {
    fn step_1(inst : &'a CSIDHInstance<K>, gen_param : i32, secret_wanted : usize) -> Bob<'a>{
        let key = inst.sample_keys(gen_param);

        Bob{
            pk: key.0,
            sk: key.1,

            k: secret_wanted,

            inst
        }
    }

    fn step_2(&self, pk_alices : [PublicKey; 2]) -> PublicKey{
        self.inst.class_group_action(pk_alices[self.k].clone(), self.sk.clone())
    }

    fn step_3_retrieve(&self, enc_mess : [u64; 2]) -> u64{
        enc_mess[self.k] ^ hash_k(self.pk.clone())
    }

}


pub fn oblivious_transfert_demo(gen_param : i32, s : [u64; 2]){

    let inst : CSIDHInstance<K> = CSIDHInstance::new(11,
        Integer::from_str_radix("14841476269619", 10).unwrap(),
        vec![Integer::from(3),Integer::from(5),Integer::from(7),
                Integer::from(11),Integer::from(13),Integer::from(17),
                Integer::from(19),Integer::from(23),Integer::from(29),
                Integer::from(31),Integer::from(37)]
    );
    inst.check_well_defined();

    println!("Step 1 -- Sampling keys");
    let alice = Alice::step_1(&inst, gen_param);

    let k = rand::random::<bool>() as usize;
    let bob = Bob::step_1(&inst, gen_param, k);

    println!("Alice's public keys:\n{}\n{}\n", alice.pk0, alice.pk1);
    println!("Bob's public key:\n{}", bob.pk);
    println!("Secret wanted by Bob: {}\n", k);

    println!("Step 2 -- Sharing secret");

    let shared_secret = bob.step_2([alice.pk0.clone(), alice.pk1.clone()]);
    println!("Shared secret:\n{}\n", shared_secret);

    let enc_mess = alice.step_2(s, shared_secret);
    println!("Encrypted messages sent by Alice:\n{}\n{}\n", enc_mess[0], enc_mess[1]);

    println!("Step 3 -- retrieving the message");
    let retrieved = bob.step_3_retrieve(enc_mess);

    println!("Retrieved message: {}\nMessages at the beggining: {}, {}", retrieved, s[0], s[1]);

}