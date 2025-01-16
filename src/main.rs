use rand::seq::SliceRandom;
use rand::Rng;


// نئوکلیوتیدهای DNA
const NUCLEOTIDES: [char; 4] = ['A', 'C', 'G', 'T'];

// تابع مقایسه
fn compare(p: char, sm: char) -> u32 {
    if p == sm { 1 } else { 0 }
}                                                                                          

// تابع جهش
fn mutate(motif: &mut Vec<char>) {
    let mut rng = rand::thread_rng();
    let index = rng.gen_range(0..motif.len());
    motif[index] = *NUCLEOTIDES.choose(&mut rng).unwrap();
}

// تابع اضافه‌کردن نوکلئوتید به موتیف
fn add_nucleotide(motif: &mut Vec<char>) {
    let mut rng = rand::thread_rng();
    motif.push(*NUCLEOTIDES.choose(&mut rng).unwrap());
}

// تابع حذف نوکلئوتید
fn delete_nucleotide(motif: &mut Vec<char>) {
    motif.pop();
}

// تابع محاسبه امتیاز فیتنس
fn fitness_score(motif: &Vec<char>, sequences: &Vec<String>) -> u32 {
    let mut score : u32 = 0;
    let motif_len = motif.len();
    for sequence in sequences {
        let mut max_match: u32 = 0;
        for i in 0..(sequence.len() - motif_len + 1) {
            let subseq: Vec<char> = sequence.chars().skip(i).take(motif_len).collect();
            let match_score: u32 = (0..motif_len).map(|j| compare(motif[j], subseq[j])).sum();
            max_match = max_match.max(match_score);
        }
        score += max_match;
    }
    score
}

fn sd_calculation(fin_store: &mut Vec<(Vec<char>, u32)>, max_length: usize,) -> Vec<u32> {
    let mut sd_score : Vec<u32> = vec![0; max_length-1];

    for i in 4..(max_length-1){
        sd_score[i] = fin_store[i+1].1 - (2*fin_store[i].1) + fin_store[i-1].1;
    }

    sd_score
}

// الگوریتم اصلی AMDILM
fn amdilm(sequences: Vec<String>, max_length: usize, max_iterations: usize) -> Vec<(Vec<char>, u32)> {
    let mut current_length = 3;
    let mut population: Vec<Vec<char>> = (0..64).map(|_| {
        let mut rng = rand::thread_rng();
        (0..current_length).map(|_| *NUCLEOTIDES.choose(&mut rng).unwrap()).collect()
    }).collect();

    let mut scores : [u32; 64] = [0; 64];
    for (i,individual) in population.iter_mut().enumerate() {
        scores[i] = fitness_score(individual, &sequences);
    }

    let mut fin_store : Vec<(Vec<char>, u32)> = vec![(vec![], 0); max_length];
    while current_length < max_length-1 {
        for (i, individual) in population.iter_mut().enumerate() {
            // عملیات اضافه کردن
            let mut extended_motif = individual.clone();
            add_nucleotide(&mut extended_motif);
            let mut score_extended = fitness_score(&extended_motif, &sequences);

            let mut mutate_extended = extended_motif.clone();
            // انجام جهش و حذف
            for _ in 0..max_iterations {
                delete_nucleotide(&mut mutate_extended);
                if current_length > 3 {
                    mutate(&mut mutate_extended);
                }
                
                add_nucleotide(&mut mutate_extended);

                let score_mutated_extended = fitness_score(&mutate_extended, &sequences);

                // نگه‌داشتن موتیفی با امتیاز بهتر
                if score_mutated_extended > score_extended {
                    extended_motif = mutate_extended.clone();
                    score_extended = score_mutated_extended;
                }
            }
            *individual = extended_motif;
            scores[i] = score_extended; 
        }

        let (maxstore_in_len, mut index_maxstore_in_len): (u32, usize) = (0,0);
        for (i, st) in scores.iter().enumerate(){
        if maxstore_in_len < *st {
            index_maxstore_in_len = i;
        }

        fin_store[current_length+1] = (population[index_maxstore_in_len].clone(), maxstore_in_len);
    }
        
        current_length += 1;    
    }

    // پیدا کردن بهترین موتیف
    

    fin_store
}

fn main() {
    let sequences = vec![
        "ACGTACGTACGT".to_string(),
        "TGCATGCATGCA".to_string(),
        "ACGTACGTACGT".to_string(),
    ];

    let mut fin_store: Vec<(Vec<char>, u32)> = amdilm(sequences, 10, 50);

    let sd_score : Vec<u32> = sd_calculation(&mut fin_store, 10);
    
    let min = sd_score.iter().copied().min().unwrap_or(0) as usize;
    
    let optimal_motif_str: String = fin_store[min].0.iter().collect();
    println!("Optimal Motif: {}, Score: {}", optimal_motif_str, min);
}