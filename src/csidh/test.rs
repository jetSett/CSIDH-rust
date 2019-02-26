use super::*;

use elliptic_curve_algorithms::finite_fields::*;

elliptic_curve_algorithms::declare_finite_field!(K, Integer, 
            Integer::from_str_radix("37118532150319619", 10).unwrap(), 
            _m);


type PublicKey = K;
type SecretKey = Vec<i32>;


fn get_global_instance() -> CSIDHInstance<K>{
    CSIDHInstance::new(13,
            Integer::from_str_radix("37118532150319619", 10).unwrap(),
            vec![Integer::from(3),Integer::from(5),Integer::from(7),
                    Integer::from(11),Integer::from(13),Integer::from(17),
                    Integer::from(19),Integer::from(23),Integer::from(29),
                    Integer::from(31),Integer::from(37),Integer::from(41),
                    Integer::from(61)]
    )
}

fn sample_montgomery() -> EllipticCurve<K>{

    let inst = get_global_instance();

    let mut a = K::new(Integer::sample_uniform(&Integer::from(0), &Integer::from(inst.p.clone())));
    let mut ell = EllipticCurve::new_montgomery(a.clone());
    while ! (inst.is_supersingular(&ell) && 
            ell.discriminant() != K::from_int(0)) {
        a = K::new(Integer::sample_uniform(&Integer::from(0), &Integer::from(inst.p.clone())));
        ell = EllipticCurve::new_montgomery(a);
    }
    ell
}

#[test]
fn everything_ok() {
    let inst = get_global_instance();

    inst.check_well_defined();
}

#[test]
fn known_supersingular_ec(){

    let inst = get_global_instance();

    // Note P = 3 mod 4, 2 mod 3
    assert!(inst.is_supersingular(&EllipticCurve::new_reduced_weierstrass(K::from_int(1), K::from_int(0)))); // supersingular iff P = 3 mod 4
    assert!(inst.is_supersingular(&EllipticCurve::new_reduced_weierstrass(K::from_int(0), K::from_int(1)))); // supersingular iff P = 2 mod 3
}

#[test]
fn class_group_action_keep_valid() {
    let inst = get_global_instance();

    for _i in 0..10{
        let mut pk = PublicKey::from_int(0);
    
        let mut action : SecretKey = vec![0, 0, 0, 0, 0, 0];
        for j in 0..6{
            action[j] = rand::random::<i32>()%3;
        }
        assert!(inst.verify_public_key(inst.class_group_action(pk, action)));
    }

}

#[test]
fn class_group_action_commute(){
    let inst = get_global_instance();
    let N_PRIMES = inst.n_primes;

    for _i in 0..3{
        let mut pk = K::from_int(0);
    
        let mut action1 : SecretKey = vec!();
        let mut action2 : SecretKey = vec!();

        for _j in 0..N_PRIMES{
            action1.push(rand::random::<i32>()%2);
            action2.push(rand::random::<i32>()%2);
        }


        let pk1 = inst.class_group_action(pk.clone(), action1.clone());
        let pk12 = inst.class_group_action(pk1.clone(), action2.clone());


        let pk2 = inst.class_group_action(pk, action2.clone());
        let pk21 = inst.class_group_action(pk2, action1.clone());

        assert_eq!(pk12, pk21);
    }
}

#[test]
fn naive_class_group_action_commute(){
    let inst = get_global_instance();
    let N_PRIMES = inst.n_primes;

    for _i in 0..10{
        let mut pk = K::from_int(0);
    
        let mut action1 : SecretKey = vec!();
        let mut action2 : SecretKey = vec!();

        for _j in 0..N_PRIMES{
            action1.push(rand::random::<i32>()%2);
            action2.push(rand::random::<i32>()%2);
        }


        let pk1 = inst.naive_class_group_action(pk.clone(), action1.clone());
        let pk12 = inst.naive_class_group_action(pk1.clone(), action2.clone());


        let pk2 = inst.naive_class_group_action(pk, action2.clone());
        let pk21 = inst.naive_class_group_action(pk2, action1.clone());

        assert_eq!(pk12, pk21);
    }
}

#[test]
fn isogeny_kernel(){
    let inst = get_global_instance();
    let N_PRIMES = inst.n_primes;

    let ell = EllipticCurve::new_montgomery(K::from_int(0));
    let mut k = Integer::from(4);

    for i in 0..N_PRIMES{
        k *= inst.l[i as usize].clone();
    }

    for _i in 0..10{
        let j = rand::random::<usize>()%(N_PRIMES as usize);
        let mut p = EllipticCurve::unsigne_point(ell.sample_point());

        p = ell.scalar_mult_unsigned(k.clone()/inst.l[j].clone(), p);
        if p == UnsignedProjPoint::infinite_point() || p.x == K::from_int(0){
            continue;
        }

        let (_, q) = CSIDHInstance::isogeny(&ell, &p, p.clone(), inst.l[j].clone()).unwrap();
        assert_eq!(UnsignedProjPoint::infinite_point(), q);
    }
}
