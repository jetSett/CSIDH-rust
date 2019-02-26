#![feature(type_alias_enum_variants)]

use time;

mod csidh;

use elliptic_curve_algorithms::finite_fields::*;

use csidh::*;

elliptic_curve_algorithms::declare_finite_field!(K, Integer, 
            Integer::from_str_radix("5326738796327623094747867617954605554069371494832722337612446642054009560026576537626892113026381253624626941643949444792662881241621373288942880288065659", 10).unwrap(), 
            _m);

type PublicKey = K;
type SecretKey = Vec<i32>;

fn main() {
    let gen_param = 1;

    let mut current : u64;

    let inst : CSIDHInstance<K> = CSIDHInstance::new(74,
        Integer::from_str_radix("5326738796327623094747867617954605554069371494832722337612446642054009560026576537626892113026381253624626941643949444792662881241621373288942880288065659", 10).unwrap(),
        vec![Integer::from(3),Integer::from(5),Integer::from(7),
                Integer::from(11),Integer::from(13),Integer::from(17),
                Integer::from(19),Integer::from(23),Integer::from(29),
                Integer::from(31),Integer::from(37),Integer::from(41),
                Integer::from(43),Integer::from(47),Integer::from(53),
                Integer::from(59),Integer::from(61),Integer::from(67),
                Integer::from(71),Integer::from(73),Integer::from(79),
                Integer::from(83),Integer::from(89),Integer::from(97),
                Integer::from(101),Integer::from(103),Integer::from(107),
                Integer::from(109),Integer::from(113),Integer::from(127),
                Integer::from(131),Integer::from(137),Integer::from(139),
                Integer::from(149),Integer::from(151),Integer::from(157),
                Integer::from(163),Integer::from(167),Integer::from(173),
                Integer::from(179),Integer::from(181),Integer::from(191),
                Integer::from(193),Integer::from(197),Integer::from(199),
                Integer::from(211),Integer::from(223),Integer::from(227),
                Integer::from(229),Integer::from(233),Integer::from(239),
                Integer::from(241),Integer::from(251),Integer::from(257),
                Integer::from(263),Integer::from(269),Integer::from(271),
                Integer::from(277),Integer::from(281),Integer::from(283),
                Integer::from(293),Integer::from(307),Integer::from(311),
                Integer::from(313),Integer::from(317),Integer::from(331),
                Integer::from(337),Integer::from(347),Integer::from(349),
                Integer::from(353),Integer::from(359),Integer::from(367),
                Integer::from(373),Integer::from(587)]
    );

    inst.check_well_defined();
    println!("Everything is well defined");

    current = time::precise_time_ns();
    let (pk_a, sk_a) = inst.sample_keys(gen_param);
    let time_gen_a = (time::precise_time_ns()-current)/1000000;
    println!("Alice's public key: {} ({} ms)", pk_a, time_gen_a);

    current = time::precise_time_ns();
    let (pk_b, sk_b) = inst.sample_keys(gen_param);
    let time_gen_b = (time::precise_time_ns()-current)/1000000;
    println!("Bob's public key: {} ({} ms)", pk_b, time_gen_b);

    current = time::precise_time_ns();
    let shared_a = inst.class_group_action(pk_b, sk_a);
    let time_shared_a = (time::precise_time_ns()-current)/1000000;
    println!("Alice's shared secret: {} ({} ms)", shared_a, time_shared_a);


    current = time::precise_time_ns();
    let shared_b = inst.class_group_action(pk_a, sk_b);
    let time_shared_b = (time::precise_time_ns()-current)/1000000;
    println!("Bob's shared secret: {} ({} ms)", shared_b, time_shared_b);

    assert_eq!(shared_a, shared_b);
    println!("Shared secret match!");
}