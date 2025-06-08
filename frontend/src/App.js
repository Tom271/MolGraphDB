import React, { useState } from 'react';
import './App.css';
import MolecularVisualization from './MolecularVisualization';
import Docs from './Docs';

function App() {
  const [smiles1, setSmiles1] = useState('');
  const [smiles2, setSmiles2] = useState('');
  const [results, setResults] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const [showDocs, setShowDocs] = useState(false);

  const handleSubmit = async (e) => {
    e.preventDefault();
    
    if (!smiles1.trim() || !smiles2.trim()) {
      setError('Please enter both SMILES codes');
      return;
    }

    setLoading(true);
    setError('');
    setResults(null);
    
    try {
      const response = await fetch('http://localhost:8000/api/calculate_similarity', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          smiles1: smiles1.trim(),
          smiles2: smiles2.trim()
        }),
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'Failed to calculate similarity');
      }

      const data = await response.json();
      setResults(data);
    } catch (err) {
      setError(err.message);
    } finally {
      setLoading(false);
    }
  };

  const handleClear = () => {
    setSmiles1('');
    setSmiles2('');
    setResults(null);
    setError('');
  };

  const loadExample = () => {
    setSmiles1('CCO'); // Ethanol
    setSmiles2('CCCO'); // Propanol
  };

  return (
    <div className="App">
      <div className="container">
        <header className="header">
          <button 
            className="docs-button"
            onClick={() => setShowDocs(true)}
            title="View Documentation"
          >
            📚 Docs
          </button>
          <div className="brand-subtitle">NEarest MOlecule</div>
          <h1>🧬 NeMo</h1>
          <p>Advanced molecular similarity analysis using graph edit distance and SMILES notation</p>
        </header>

        <form onSubmit={handleSubmit} className="form">
          <div className="input-group">
            <div className="input-field">
              <label htmlFor="smiles1">Molecule 1 (SMILES):</label>
              <input
                type="text"
                id="smiles1"
                value={smiles1}
                onChange={(e) => setSmiles1(e.target.value)}
                placeholder="e.g., CCO (ethanol)"
                className="input"
              />
            </div>
            
            <div className="input-field">
              <label htmlFor="smiles2">Molecule 2 (SMILES):</label>
              <input
                type="text"
                id="smiles2"
                value={smiles2}
                onChange={(e) => setSmiles2(e.target.value)}
                placeholder="e.g., CCCO (propanol)"
                className="input"
              />
            </div>
          </div>

          <div className="button-group">
            <button
              type="submit"
              disabled={loading}
              className="btn btn-primary"
            >
              {loading ? 'Calculating...' : 'Calculate Similarity'}
            </button>
            
            <button
              type="button"
              onClick={loadExample}
              className="btn btn-secondary"
            >
              Load Example
            </button>
            
            <button
              type="button"
              onClick={handleClear}
              className="btn btn-tertiary"
            >
              Clear
            </button>
          </div>
        </form>

        {error && (
          <div className="error">
            <p>❌ {error}</p>
          </div>
        )}

        {results && (
          <div className="results">
            <h2>📊 Results</h2>
            <div className="results-grid">
              <div className="result-card">
                <h3>Similarity Score</h3>
                <div className="result-value">
                  {(results.similarity_score * 100).toFixed(1)}%
                </div>
              </div>
              
              <div className="result-card">
                <h3>Graph Edit Distance</h3>
                <div className="result-value">
                  {results.graph_edit_distance.toFixed(2)}
                </div>
              </div>
              
              <div className="result-card">
                <h3>Subgraphs Found</h3>
                <div className="result-value">
                  Mol1: {results.subgraphs1_count} | Mol2: {results.subgraphs2_count}
                </div>
              </div>
            </div>

            {/* Molecule Images */}
            {(results.molecule1_viz?.image || results.molecule2_viz?.image) && (
              <div className="molecule-images">
                <h3>🧬 Molecular Structures</h3>
                <div className="image-grid">
                  {results.molecule1_viz?.image && (
                    <div className="molecule-image-card">
                      <h4>Molecule 1</h4>
                      <img 
                        src={results.molecule1_viz.image} 
                        alt="Molecule 1 structure" 
                        className="molecule-image"
                      />
                      <p className="molecule-smiles">{results.molecule1_info.smiles}</p>
                    </div>
                  )}
                  {results.molecule2_viz?.image && (
                    <div className="molecule-image-card">
                      <h4>Molecule 2</h4>
                      <img 
                        src={results.molecule2_viz.image} 
                        alt="Molecule 2 structure" 
                        className="molecule-image"
                      />
                      <p className="molecule-smiles">{results.molecule2_info.smiles}</p>
                    </div>
                  )}
                </div>
              </div>
            )}

            <div className="molecule-info">
              <div className="molecule-card">
                <h4>Molecule 1 Details</h4>
                <p><strong>SMILES:</strong> {results.molecule1_info.smiles}</p>
                <p><strong>Formula:</strong> {results.molecule1_info.molecular_formula}</p>
                <p><strong>Atoms:</strong> {results.molecule1_info.num_atoms}</p>
                <p><strong>Bonds:</strong> {results.molecule1_info.num_bonds}</p>
                <p><strong>Subgraphs:</strong> {results.molecule1_viz?.subgraphs?.length || 0}</p>
              </div>
              
              <div className="molecule-card">
                <h4>Molecule 2 Details</h4>
                <p><strong>SMILES:</strong> {results.molecule2_info.smiles}</p>
                <p><strong>Formula:</strong> {results.molecule2_info.molecular_formula}</p>
                <p><strong>Atoms:</strong> {results.molecule2_info.num_atoms}</p>
                <p><strong>Bonds:</strong> {results.molecule2_info.num_bonds}</p>
                <p><strong>Subgraphs:</strong> {results.molecule2_viz?.subgraphs?.length || 0}</p>
              </div>
            </div>

            {/* Molecular Visualization */}
            <MolecularVisualization results={results} />
          </div>
        )}

        <footer className="footer">
          <p>💡 <strong>NeMo Tip:</strong> Use simple carbon-based molecules (e.g., CCO, CCCC, C1CCC1) for optimal molecular similarity analysis</p>
        </footer>
      </div>
      
      {showDocs && <Docs onClose={() => setShowDocs(false)} />}
    </div>
  );
}

export default App;