import React from 'react';
import './Docs.css';

const Docs = ({ onClose }) => {
  return (
    <div className="docs-overlay" onClick={onClose}>
      <div className="docs-container" onClick={(e) => e.stopPropagation()}>
        <div className="docs-header">
          <h2>🧬 NeMo Documentation</h2>
          <button className="docs-close" onClick={onClose}>×</button>
        </div>
        
        <div className="docs-content">
          <section className="docs-section">
            <h3>About NeMo</h3>
            <p>
              NeMo (NEarest MOlecule) is a molecular similarity analysis tool that compares molecules 
              using graph edit distance. It visualizes molecular subgraphs and their relationships 
              to help understand structural similarities between compounds.
            </p>
          </section>

          <section className="docs-section">
            <h3>How to Use</h3>
            <div className="docs-steps">
              <div className="step">
                <strong>1.</strong> Enter two SMILES strings representing your molecules
              </div>
              <div className="step">
                <strong>2.</strong> Click "Calculate Similarity" to analyze the molecules
              </div>
              <div className="step">
                <strong>3.</strong> View the similarity score and visualization below
              </div>
            </div>
            <p className="docs-note">
              <strong>Note:</strong> Use simple carbon-based molecules (e.g., CCO, CCCC, C1CCC1) for best results.
            </p>
          </section>

          <section className="docs-section">
            <h3>Graph Edit Distance Algorithm</h3>
            <p>
              The Maximum Common Substructure (MCS) algorithm works by:
            </p>
            <ul className="docs-list">
              <li><strong>Subgraph Generation:</strong> Breaking molecules into all possible connected subgraphs</li>
              <li><strong>Canonical Hashing:</strong> Creating unique identifiers for each subgraph structure</li>
              <li><strong>Overlap Analysis:</strong> Finding shared subgraphs between molecules using set intersection</li>
              <li><strong>Distance Calculation:</strong> Computing similarity based on Tanimoto coefficient and shared structures</li>
            </ul>
          </section>

          <section className="docs-section">
            <h3>Visualization</h3>
            <p>
              The tree visualization shows:
            </p>
            <ul className="docs-list">
              <li><span className="color-box mol1"></span> <strong>Brown:</strong> Molecule 1 subgraphs</li>
              <li><span className="color-box mol2"></span> <strong>Green:</strong> Molecule 2 subgraphs</li>
              <li><span className="color-box shared"></span> <strong>Orange:</strong> Shared subgraphs</li>
            </ul>
            <p>
              Subgraphs are organized by edge count (complex to simple, top to bottom) with 
              hierarchical connections showing parent-child relationships.
            </p>
          </section>

          <section className="docs-section">
            <h3>Technical Details</h3>
            <div className="tech-grid">
              <div className="tech-item">
                <strong>Backend:</strong> FastAPI with RDKit for molecule processing
              </div>
              <div className="tech-item">
                <strong>Frontend:</strong> React with SVG-based visualization
              </div>
              <div className="tech-item">
                <strong>Algorithm:</strong> NetworkX for graph operations and MCS calculation
              </div>
              <div className="tech-item">
                <strong>Database:</strong> SQLite for subgraph storage and caching
              </div>
            </div>
          </section>
        </div>
      </div>
    </div>
  );
};

export default Docs;