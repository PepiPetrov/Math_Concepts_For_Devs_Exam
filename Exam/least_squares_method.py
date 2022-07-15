import numpy as np

def least_sq(sample_spectrum, components):
    
    # Sample_spectrum (unknown spectrum): array of w values.
    # Components (known spectra): array of n (number of components) columns with w values.
    # This def returns an array of n values. Each value is the similarity score for the sample_spectrum and a component spectrum.

    similarity = np.dot(np.linalg.inv(np.dot(components, components.T)) , np.dot(components, sample_spectrum))
    
    return similarity