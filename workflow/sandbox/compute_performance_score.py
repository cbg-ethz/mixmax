from pathlib import Path
import logging

import yaml
import pandas as pd
from sklearn.metrics import v_measure_score
from statistics import mean

logging.basicConfig(format="{asctime} - {levelname} - {message}", style="{", datefmt="%Y-%m-%d %H:%M",level=logging.DEBUG)


def compute_score(labels, estimate):
    mismatches = 0
    total_counts = 0
    for key1, dataset2sampleID in estimate.items():
        assignments = dataset2sampleID.items()
        total_counts += len(assignments)
        for key2, sample_assignment in assignments:
            if sample_assignment != labels[str(key1)][str(key2)]:
                mismatches += 1
            
    return mismatches / total_counts
#def preprocess_output():



def load_data(path):
    with open(path / 'sample_identity.yaml', 'r') as file:
        sample_identity = yaml.safe_load(file)

    with open(path / 'sample_assignment.yaml', 'r') as file:
        sample_assignment = yaml.safe_load(file)
    
    return sample_identity, sample_assignment


def preprocess_labels(sample_identity):
    for key, value_dict in sample_identity.items():
        values = [v for k, v in value_dict.items() if k != '-1']
        for sub_key, value in value_dict.items():
            value_dict[sub_key] = int(value.replace('+', ''))
        if len(values) != len(set(values)):
            logging.error(f"Duplicate values found in sample_identity for key {key}: {values}")
            return True
    return False



def preprocess_sample_assignment(sample_assignment):
    for key, value_dict in sample_assignment.items():
        for sub_key, value in value_dict.items():
            value_dict[sub_key] = int(value.split('-')[1])


def process_simulation_run(path):
    sample_identity, sample_assignment = load_data(path)
    sample_identity_is_pathological = preprocess_labels(sample_identity)
    preprocess_sample_assignment(sample_assignment)
    
    logging.info(sample_identity_is_pathological)
    score = compute_score(sample_identity, sample_assignment)

    return sample_identity_is_pathological, score



def compute_clustering_performance_score(subfolder):
    demultiplexed_folder = subfolder / 'demultiplexed'

    v_measures = []
    for assignment_file in demultiplexed_folder.glob('*assignments.tsv'):
        df = pd.read_csv(assignment_file, sep='\t')
        
        # Extract ground truth labels and cluster assignments

        ground_truth_labels = df.columns[1:].map(lambda col: 'doublet' if '+' in col else col.split('_')[1].split('-')[1]).tolist()
        cluster_assignments = df.iloc[0, 1:].tolist()

        # Compute ARI

        v_measure = v_measure_score(ground_truth_labels, cluster_assignments)
        v_measures.append(v_measure)

    return mean(v_measures)


""" def main(input_dir):
    if not isinstance(input_dir, Path):
        input_dir = Path(input_dir)
    results = []
    for seed_folder in input_dir.glob('seed*'):
        for subfolder_0 in seed_folder.glob('0*'):
            for subfolder in subfolder_0.glob('*0'):
                sample_identity_path = subfolder / 'sample_identity.yaml'
                sample_assignment_path = subfolder / 'sample_assignment.yaml'
                
                if sample_identity_path.exists() and sample_assignment_path.exists():
                    sample_identity_is_pathological, score = process_simulation_run(subfolder)

                    mean_v_score = compute_clustering_performance_score(subfolder)

                    seed_number = int(seed_folder.name.split('_')[1])
                    
                    results.append({
                        'seed': seed_number,
                        'doublet_rate': float(subfolder_0.name),
                        'cell_count': int(subfolder.name),
                        'score': score,
                        'sample_identity_is_pathological': sample_identity_is_pathological,
                        'v_score': mean_v_score
                    })

                    # Periodically save results to a file to avoid memory burden
                    if len(results) % 100 == 0:
                        temp_df = pd.DataFrame(results)
                        temp_df.to_csv('more_intermediate_results.csv', mode='a', header=False, index=False)
                        results.clear() """

def main(input_dir):
    if not isinstance(input_dir, Path):
        input_dir = Path(input_dir)
    results = []

    for seed_folder in input_dir.glob('seed*'):
        subfolder = seed_folder / "robust~True" / "pool_size~4" / "doublet_rate~0.05" / "cell_count~1000"
        sample_assignment_path = subfolder /  "sample_assignment.yaml" 
        sample_identity_path = subfolder / "sample_identity.yaml" 
        
        logging.info(f"Processing simulation run for {subfolder}")

        if sample_identity_path.exists() and sample_assignment_path.exists():
            logging.info("Subfolder exists!")
            sample_identity_is_pathological, score = process_simulation_run(subfolder)

            mean_v_score = compute_clustering_performance_score(subfolder)

            seed_number = int(seed_folder.name.split('_')[1])
                    
            results.append({
                'seed': seed_number,
                'doublet_rate': 0.05,
                'cell_count': 1000,
                'score': score,
                'sample_identity_is_pathological': sample_identity_is_pathological,
                'v_score': mean_v_score
            })

                    # Periodically save results to a file to avoid memory burden
            if len(results) % 100 == 0:
                temp_df = pd.DataFrame(results)
                temp_df.to_csv('/cluster/work/bewi/members/jgawron/projects/Demultiplexing/Demultiplexing_simulations/workflow/sandbox/results_for_robust.csv', mode='a', header=False, index=False)
                results.clear()


    # Save any remaining results
    if results:
        df = pd.DataFrame(results)
        df.to_csv('/cluster/work/bewi/members/jgawron/projects/Demultiplexing/Demultiplexing_simulations/workflow/sandbox/results_for_robust.csv', mode='a', header=not Path('/cluster/work/bewi/members/jgawron/projects/Demultiplexing/Demultiplexing_simulations/workflow/sandbox/results_for_robust.csv').exists(), index=False)


if __name__ == '__main__':
    input_dir = Path('/cluster/work/bewi/members/jgawron/projects/Demultiplexing/AML_data/output')# = args.input
    main(input_dir)
