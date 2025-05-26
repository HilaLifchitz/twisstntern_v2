import numpy as np
import pandas as pd
import twisstntern

def generate_sample_data(n_samples=1000):
    """Generate sample data for testing."""
    np.random.seed(42)
    t1 = np.random.normal(0.3, 0.1, n_samples)
    t2 = np.random.normal(0.5, 0.1, n_samples)
    t3 = np.random.normal(0.7, 0.1, n_samples)
    
    # Ensure values are between 0 and 1
    t1 = np.clip(t1, 0, 1)
    t2 = np.clip(t2, 0, 1)
    t3 = np.clip(t3, 0, 1)
    
    return pd.DataFrame({
        'T1': t1,
        'T2': t2,
        'T3': t3
    })

def main():
    # Generate sample data
    print("Generating sample data...")
    data = generate_sample_data()
    
    # Save sample data to CSV
    data.to_csv('sample_data.csv', index=False)
    print("Sample data saved to 'sample_data.csv'")
    
    # Run analysis with default parameters
    print("\nRunning analysis...")
    twisstntern.run_analysis('sample_data.csv', granuality='superfine')
    
    print("\nTest completed successfully!")

if __name__ == "__main__":
    main() 