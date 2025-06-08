import React, { useEffect, useRef } from 'react';
import './MolecularVisualization.css';

const MolecularVisualization = ({ results }) => {
  const svgRef = useRef(null);

  // Calculate dynamic height based on content
  const calculateHeight = (results) => {
    if (!results?.molecule1_viz || !results?.molecule2_viz) return 800;
    
    // Group by edge count to estimate levels
    const allEdgeCounts = new Set();
    results.molecule1_viz.subgraphs.forEach(sg => allEdgeCounts.add(sg.edge_count));
    results.molecule2_viz.subgraphs.forEach(sg => allEdgeCounts.add(sg.edge_count));
    
    const numLevels = allEdgeCounts.size;
    const baseHeight = 400; // For molecules at top
    const levelHeight = 140;
    const paddingHeight = 150;
    
    return Math.max(800, baseHeight + (numLevels * levelHeight) + paddingHeight);
  };

  // Configuration for tree layout
  const config = {
    width: 1400,
    height: calculateHeight(results),
    margin: { top: 50, right: 50, bottom: 50, left: 50 },
    moleculeSize: { width: 300, height: 250 },
    subgraphSize: { width: 100, height: 80 },
    levelHeight: 140,
    nodeSpacing: 120,
  };

  useEffect(() => {
    if (!results?.molecule1_viz || !results?.molecule2_viz) return;

    const svg = svgRef.current;
    if (!svg) return;

    // Clear previous content
    svg.innerHTML = '';

    // Create SVG namespace
    const svgNS = 'http://www.w3.org/2000/svg';

    // Set SVG dimensions
    svg.setAttribute('width', config.width);
    svg.setAttribute('height', config.height);
    svg.setAttribute('viewBox', `0 0 ${config.width} ${config.height}`);

    // Create tree layout data
    const treeData = createTreeLayout(results);

    // Draw connections first (so they appear behind nodes)
    drawConnections(svg, treeData, svgNS);

    // Draw nodes
    drawNodes(svg, treeData, svgNS);

    // Add legend
    drawLegend(svg, svgNS);

  }, [results]);

  const createTreeLayout = (results) => {
    const mol1 = results.molecule1_viz;
    const mol2 = results.molecule2_viz;
    const sharedIds = new Set(results.shared_subgraphs);

    // Group subgraphs by edge count for tree levels
    const groupByEdgeCount = (subgraphs) => {
      const groups = {};
      subgraphs.forEach(sg => {
        if (!groups[sg.edge_count]) {
          groups[sg.edge_count] = [];
        }
        groups[sg.edge_count].push(sg);
      });
      return groups;
    };

    const mol1Groups = groupByEdgeCount(mol1.subgraphs);
    const mol2Groups = groupByEdgeCount(mol2.subgraphs);

    // Get all edge counts and sort them in descending order
    const allEdgeCounts = new Set([
      ...Object.keys(mol1Groups).map(Number),
      ...Object.keys(mol2Groups).map(Number)
    ]);
    const sortedEdgeCounts = Array.from(allEdgeCounts).sort((a, b) => b - a);

    // Calculate layout positions
    const layout = {
      molecules: [],
      subgraphs: [],
      connections: []
    };

    // Position main molecules at the top
    const moleculeY = config.margin.top;
    const mol1X = config.width * 0.25;
    const mol2X = config.width * 0.75;

    layout.molecules.push({
      type: 'molecule',
      data: mol1,
      x: mol1X,
      y: moleculeY,
      id: 'mol1'
    });

    layout.molecules.push({
      type: 'molecule',
      data: mol2,
      x: mol2X,
      y: moleculeY,
      id: 'mol2'
    });

    // Position subgraphs in tree levels
    sortedEdgeCounts.forEach((edgeCount, levelIndex) => {
      const levelY = moleculeY + config.moleculeSize.height + (levelIndex + 1) * config.levelHeight;
      
      const mol1Subgraphs = mol1Groups[edgeCount] || [];
      const mol2Subgraphs = mol2Groups[edgeCount] || [];
      
      // Create a merged list for shared subgraphs
      const allSubgraphsAtLevel = [];
      const processedShared = new Set();

      // Add mol1 subgraphs
      mol1Subgraphs.forEach(sg => {
        if (sharedIds.has(sg.id) && !processedShared.has(sg.id)) {
          // This is a shared subgraph - position it in the center
          allSubgraphsAtLevel.push({
            ...sg,
            source: 'shared',
            parentMol: 'both'
          });
          processedShared.add(sg.id);
        } else if (!sharedIds.has(sg.id)) {
          // This is unique to mol1
          allSubgraphsAtLevel.push({
            ...sg,
            source: 'mol1',
            parentMol: 'mol1'
          });
        }
      });

      // Add mol2 subgraphs (only unique ones, shared ones are already added)
      mol2Subgraphs.forEach(sg => {
        if (!sharedIds.has(sg.id)) {
          allSubgraphsAtLevel.push({
            ...sg,
            source: 'mol2',
            parentMol: 'mol2'
          });
        }
      });

      // Calculate positions for this level
      const totalWidth = config.width - 2 * config.margin.left;
      const availableWidth = totalWidth - config.subgraphSize.width;
      const spacing = allSubgraphsAtLevel.length > 1 
        ? availableWidth / (allSubgraphsAtLevel.length - 1) 
        : 0;

      allSubgraphsAtLevel.forEach((sg, index) => {
        let x;
        if (allSubgraphsAtLevel.length === 1) {
          x = config.width / 2;
        } else {
          x = config.margin.left + config.subgraphSize.width / 2 + index * spacing;
        }

        const node = {
          type: 'subgraph',
          data: sg,
          x: x,
          y: levelY,
          id: sg.id,
          level: levelIndex,
          edgeCount: edgeCount,
          source: sg.source
        };

        layout.subgraphs.push(node);

        // Store node for later connection processing
        // We'll add connections after all nodes are positioned
      });
    });

    // Add hierarchical connections after all nodes are positioned
    addHierarchicalConnections(layout, mol1X, mol2X, moleculeY, config);

    return layout;
  };

  const addHierarchicalConnections = (layout, mol1X, mol2X, moleculeY, config) => {
    // First, connect ALL subgraphs to their parent molecules
    layout.subgraphs.forEach(node => {
      if (node.source === 'mol1') {
        // Connect mol1 subgraphs to molecule 1
        layout.connections.push({
          from: { x: mol1X, y: moleculeY + config.moleculeSize.height },
          to: { x: node.x, y: node.y },
          type: 'mol1'
        });
      } else if (node.source === 'mol2') {
        // Connect mol2 subgraphs to molecule 2
        layout.connections.push({
          from: { x: mol2X, y: moleculeY + config.moleculeSize.height },
          to: { x: node.x, y: node.y },
          type: 'mol2'
        });
      } else if (node.source === 'shared') {
        // Connect shared subgraphs to both molecules
        layout.connections.push({
          from: { x: mol1X, y: moleculeY + config.moleculeSize.height },
          to: { x: node.x, y: node.y },
          type: 'shared'
        });
        layout.connections.push({
          from: { x: mol2X, y: moleculeY + config.moleculeSize.height },
          to: { x: node.x, y: node.y },
          type: 'shared'
        });
      }
    });

    // Then, add hierarchical connections between subgraphs
    // Group subgraphs by edge count and source
    const subgraphsByEdgeCount = {};
    layout.subgraphs.forEach(node => {
      const key = `${node.edgeCount}_${node.source}`;
      if (!subgraphsByEdgeCount[key]) {
        subgraphsByEdgeCount[key] = [];
      }
      subgraphsByEdgeCount[key].push(node);
    });

    // Sort edge counts in descending order
    const edgeCounts = [...new Set(layout.subgraphs.map(n => n.edgeCount))].sort((a, b) => b - a);

    // Connect subgraphs in hierarchical manner (parent to child)
    edgeCounts.forEach((edgeCount, levelIndex) => {
      if (levelIndex < edgeCounts.length - 1) {
        // Connect to child subgraphs (one level down with one fewer edge)
        const childEdgeCount = edgeCounts[levelIndex + 1];
        
        ['mol1', 'mol2', 'shared'].forEach(source => {
          const parentNodes = subgraphsByEdgeCount[`${edgeCount}_${source}`] || [];
          const childNodes = subgraphsByEdgeCount[`${childEdgeCount}_${source}`] || [];
          
          // For shared subgraphs, also consider connections to mol1 and mol2 children
          let allChildNodes = [...childNodes];
          if (source === 'shared') {
            allChildNodes = [
              ...childNodes,
              ...(subgraphsByEdgeCount[`${childEdgeCount}_mol1`] || []),
              ...(subgraphsByEdgeCount[`${childEdgeCount}_mol2`] || [])
            ];
          }

          parentNodes.forEach(parentNode => {
            allChildNodes.forEach(childNode => {
              // Check if child could be contained in parent
              const score = calculateSubgraphSimilarity(childNode, parentNode);
              if (score > 0) {
                const connectionType = source === 'shared' ? 'shared' : source;
                layout.connections.push({
                  from: { x: parentNode.x, y: parentNode.y + config.subgraphSize.height },
                  to: { x: childNode.x, y: childNode.y },
                  type: connectionType
                });
              }
            });
          });
        });
      }
    });
  };

  const calculateSubgraphSimilarity = (childNode, parentNode) => {
    // Check if childNode's subgraph could be contained in parentNode's subgraph
    const childNodes = new Set(childNode.data.nodes);
    const parentNodes = new Set(parentNode.data.nodes);
    const childEdges = new Set(childNode.data.edges.map(e => `${Math.min(e[0], e[1])}-${Math.max(e[0], e[1])}`));
    const parentEdges = new Set(parentNode.data.edges.map(e => `${Math.min(e[0], e[1])}-${Math.max(e[0], e[1])}`));

    // Child must be smaller than parent
    if (childNodes.size >= parentNodes.size) return 0;

    // All child nodes should be in parent
    const nodeContainment = [...childNodes].every(node => parentNodes.has(node));
    if (!nodeContainment) return 0;

    // All child edges should be in parent
    const edgeContainment = [...childEdges].every(edge => parentEdges.has(edge));
    if (!edgeContainment) return 0;

    // Calculate containment score (higher is better)
    const nodeRatio = childNodes.size / parentNodes.size;
    const edgeRatio = childEdges.size / Math.max(parentEdges.size, 1);
    
    // Prefer parents that are "closest" in size (not too much bigger)
    const sizeScore = 1.0 - Math.abs(parentNodes.size - childNodes.size - 1) / parentNodes.size;
    
    return nodeRatio * 0.4 + edgeRatio * 0.4 + sizeScore * 0.2;
  };

  const drawConnections = (svg, layout, svgNS) => {
    layout.connections.forEach(conn => {
      const line = document.createElementNS(svgNS, 'line');
      line.setAttribute('x1', conn.from.x);
      line.setAttribute('y1', conn.from.y);
      line.setAttribute('x2', conn.to.x);
      line.setAttribute('y2', conn.to.y);
      
      let strokeColor = '#cccccc';
      let strokeWidth = '2';
      
      if (conn.type === 'shared') {
        strokeColor = '#e74c3c';
        strokeWidth = '3';
      } else if (conn.type === 'mol1') {
        strokeColor = '#3498db';
      } else if (conn.type === 'mol2') {
        strokeColor = '#2ecc71';
      }
      
      line.setAttribute('stroke', strokeColor);
      line.setAttribute('stroke-width', strokeWidth);
      line.setAttribute('opacity', '0.7');
      
      svg.appendChild(line);
    });
  };

  const drawNodes = (svg, layout, svgNS) => {
    // Draw molecules
    layout.molecules.forEach(mol => {
      const group = document.createElementNS(svgNS, 'g');
      group.setAttribute('transform', `translate(${mol.x - config.moleculeSize.width/2}, ${mol.y})`);
      
      // Background rectangle
      const rect = document.createElementNS(svgNS, 'rect');
      rect.setAttribute('width', config.moleculeSize.width);
      rect.setAttribute('height', config.moleculeSize.height);
      rect.setAttribute('rx', '8');
      rect.setAttribute('fill', '#f8f9fa');
      rect.setAttribute('stroke', mol.id === 'mol1' ? '#3498db' : '#2ecc71');
      rect.setAttribute('stroke-width', '3');
      group.appendChild(rect);
      
      // Molecule image
      if (mol.data.image) {
        const image = document.createElementNS(svgNS, 'image');
        image.setAttributeNS('http://www.w3.org/1999/xlink', 'href', mol.data.image);
        image.setAttribute('x', '10');
        image.setAttribute('y', '10');
        image.setAttribute('width', '200');
        image.setAttribute('height', '150');
        group.appendChild(image);
      }
      
      // Label
      const text = document.createElementNS(svgNS, 'text');
      text.setAttribute('x', config.moleculeSize.width / 2);
      text.setAttribute('y', config.moleculeSize.height - 20);
      text.setAttribute('text-anchor', 'middle');
      text.setAttribute('font-family', 'Arial, sans-serif');
      text.setAttribute('font-size', '14');
      text.setAttribute('font-weight', 'bold');
      text.textContent = `${mol.data.smiles} (${mol.data.subgraphs.length} subgraphs)`;
      group.appendChild(text);
      
      svg.appendChild(group);
    });

    // Draw subgraphs
    layout.subgraphs.forEach(sg => {
      const group = document.createElementNS(svgNS, 'g');
      group.setAttribute('transform', `translate(${sg.x - config.subgraphSize.width/2}, ${sg.y})`);
      
      // Background rectangle
      const rect = document.createElementNS(svgNS, 'rect');
      rect.setAttribute('width', config.subgraphSize.width);
      rect.setAttribute('height', config.subgraphSize.height);
      rect.setAttribute('rx', '4');
      
      let fillColor = '#ffffff';
      let strokeColor = '#cccccc';
      let strokeWidth = '2';
      
      if (sg.source === 'shared') {
        fillColor = '#fee';
        strokeColor = '#e74c3c';
        strokeWidth = '3';
      } else if (sg.source === 'mol1') {
        fillColor = '#eff8ff';
        strokeColor = '#3498db';
      } else if (sg.source === 'mol2') {
        fillColor = '#eefff0';
        strokeColor = '#2ecc71';
      }
      
      rect.setAttribute('fill', fillColor);
      rect.setAttribute('stroke', strokeColor);
      rect.setAttribute('stroke-width', strokeWidth);
      group.appendChild(rect);
      
      // Subgraph image
      if (sg.data.image) {
        const image = document.createElementNS(svgNS, 'image');
        image.setAttributeNS('http://www.w3.org/1999/xlink', 'href', sg.data.image);
        image.setAttribute('x', '5');
        image.setAttribute('y', '5');
        image.setAttribute('width', '60');
        image.setAttribute('height', '45');
        group.appendChild(image);
      }
      
      // Edge count label
      const edgeText = document.createElementNS(svgNS, 'text');
      edgeText.setAttribute('x', config.subgraphSize.width / 2);
      edgeText.setAttribute('y', config.subgraphSize.height - 5);
      edgeText.setAttribute('text-anchor', 'middle');
      edgeText.setAttribute('font-family', 'Arial, sans-serif');
      edgeText.setAttribute('font-size', '10');
      edgeText.setAttribute('font-weight', 'bold');
      edgeText.textContent = `${sg.edgeCount} edges`;
      group.appendChild(edgeText);
      
      svg.appendChild(group);
    });
  };

  const drawLegend = (svg, svgNS) => {
    const legendGroup = document.createElementNS(svgNS, 'g');
    legendGroup.setAttribute('transform', `translate(${config.width - 200}, 20)`);
    
    // Legend background
    const legendRect = document.createElementNS(svgNS, 'rect');
    legendRect.setAttribute('width', '200');
    legendRect.setAttribute('height', '130');
    legendRect.setAttribute('fill', '#f8f9fa');
    legendRect.setAttribute('stroke', '#dee2e6');
    legendRect.setAttribute('stroke-width', '1');
    legendRect.setAttribute('rx', '4');
    legendGroup.appendChild(legendRect);
    
    // Legend title
    const title = document.createElementNS(svgNS, 'text');
    title.setAttribute('x', '10');
    title.setAttribute('y', '20');
    title.setAttribute('font-family', 'Arial, sans-serif');
    title.setAttribute('font-size', '12');
    title.setAttribute('font-weight', 'bold');
    title.textContent = 'Legend';
    legendGroup.appendChild(title);
    
    // Legend items
    const legendItems = [
      { color: '#3498db', text: 'Molecule 1 subgraphs', y: 40 },
      { color: '#2ecc71', text: 'Molecule 2 subgraphs', y: 55 },
      { color: '#e74c3c', text: 'Shared subgraphs', y: 70 },
      { color: '#666666', text: 'Hierarchical connections', y: 85 },
      { color: '#666666', text: '(parent → child subgraphs)', y: 100 }
    ];
    
    legendItems.forEach(item => {
      const circle = document.createElementNS(svgNS, 'circle');
      circle.setAttribute('cx', '20');
      circle.setAttribute('cy', item.y - 3);
      circle.setAttribute('r', '4');
      circle.setAttribute('fill', item.color);
      legendGroup.appendChild(circle);
      
      const text = document.createElementNS(svgNS, 'text');
      text.setAttribute('x', '30');
      text.setAttribute('y', item.y);
      text.setAttribute('font-family', 'Arial, sans-serif');
      text.setAttribute('font-size', '10');
      text.textContent = item.text;
      legendGroup.appendChild(text);
    });
    
    svg.appendChild(legendGroup);
  };

  if (!results?.molecule1_viz || !results?.molecule2_viz) {
    return (
      <div className="visualization-placeholder">
        <p>No visualization data available. Run a calculation to see the molecular tree visualization.</p>
      </div>
    );
  }

  return (
    <div className="molecular-visualization">
      <h3>🌳 Molecular Subgraph Tree Visualization</h3>
      <p className="visualization-description">
        This hierarchical tree shows molecular decomposition from complex (top) to simple subgraphs (bottom). 
        Connections only link subgraphs that differ by exactly one edge, creating proper parent-child relationships.
        Shared subgraphs are highlighted in red.
      </p>
      <div className="svg-container">
        <svg ref={svgRef} className="visualization-svg"></svg>
      </div>
      <div className="visualization-stats">
        <div className="stat">
          <strong>Molecule 1:</strong> {results.molecule1_viz.subgraphs.length} subgraphs
        </div>
        <div className="stat">
          <strong>Molecule 2:</strong> {results.molecule2_viz.subgraphs.length} subgraphs
        </div>
        <div className="stat">
          <strong>Shared:</strong> {results.shared_subgraphs.length} subgraphs
        </div>
      </div>
    </div>
  );
};

export default MolecularVisualization;