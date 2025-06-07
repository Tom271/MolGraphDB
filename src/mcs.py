#!/usr/bin/env python3
"""Molecular Graph Edit Distance Calculator
Inspired by NextMove SmallWorld approach for subgraph-based similarity

This project calculates graph edit distance between carbon-only molecules
using a database of molecular subgraphs for efficient computation.
"""

import json
import logging
import sqlite3
from dataclasses import dataclass
from itertools import combinations
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import networkx as nx
from pydantic import BaseModel, Field, validator
from rdkit import Chem
from rdkit.Chem import Draw

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class MoleculeInput(BaseModel):
    """Pydantic model for validating molecular input"""

    smiles: str = Field(..., description="SMILES string of the molecule")
    name: Optional[str] = Field(None, description="Optional name for the molecule")

    @validator("smiles")
    def validate_smiles(cls, v):
        """Validate SMILES string and ensure carbon-only with single bonds"""
        mol = Chem.MolFromSmiles(v)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {v}")

        # Check atom count
        if mol.GetNumAtoms() > 10:
            raise ValueError(f"Molecule has {mol.GetNumAtoms()} atoms, maximum is 10")

        # Check only carbon atoms
        for atom in mol.GetAtoms():
            if atom.GetSymbol() != "C":
                raise ValueError(f"Non-carbon atom found: {atom.GetSymbol()}")

        # Check only single bonds
        for bond in mol.GetBonds():
            if bond.GetBondType() != Chem.BondType.SINGLE:
                raise ValueError(f"Non-single bond found: {bond.GetBondType()}")

        return v


@dataclass
class SubgraphEntry:
    """Represents a subgraph entry in the database"""

    subgraph_hash: str
    size: int
    adjacency_matrix: List[List[int]]
    parent_molecule: str
    frequency: int = 1


class SubgraphDatabase:
    """Database for storing and querying molecular subgraphs"""

    def __init__(self, db_path: str = "subgraphs.db"):
        self.db_path = db_path
        self.conn = None
        self._init_database()

    def _init_database(self):
        """Initialise the SQLite database"""
        self.conn = sqlite3.connect(self.db_path)
        cursor = self.conn.cursor()

        cursor.execute("""
            CREATE TABLE IF NOT EXISTS subgraphs (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                subgraph_hash TEXT UNIQUE,
                size INTEGER,
                adjacency_matrix TEXT,
                parent_molecules TEXT,
                frequency INTEGER DEFAULT 1
            )
        """)

        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_hash ON subgraphs(subgraph_hash)
        """)
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_size ON subgraphs(size)
        """)

        self.conn.commit()

    def add_subgraph(self, entry: SubgraphEntry):
        """Add or update a subgraph entry"""
        cursor = self.conn.cursor()

        # Check if subgraph already exists
        cursor.execute(
            "SELECT id, parent_molecules, frequency FROM subgraphs WHERE subgraph_hash = ?",
            (entry.subgraph_hash,),
        )

        result = cursor.fetchone()
        if result:
            # Update existing entry
            existing_id, existing_parents, existing_freq = result
            parents_list = json.loads(existing_parents)
            if entry.parent_molecule not in parents_list:
                parents_list.append(entry.parent_molecule)

            cursor.execute(
                """
                UPDATE subgraphs 
                SET parent_molecules = ?, frequency = ?
                WHERE id = ?
            """,
                (json.dumps(parents_list), existing_freq + 1, existing_id),
            )
        else:
            # Insert new entry
            cursor.execute(
                """
                INSERT INTO subgraphs 
                (subgraph_hash, size, adjacency_matrix, parent_molecules, frequency)
                VALUES (?, ?, ?, ?, ?)
            """,
                (
                    entry.subgraph_hash,
                    entry.size,
                    json.dumps(entry.adjacency_matrix),
                    json.dumps([entry.parent_molecule]),
                    entry.frequency,
                ),
            )

        self.conn.commit()

    def get_subgraphs_by_size(self, size: int) -> List[Dict]:
        """Retrieve all subgraphs of a given size"""
        cursor = self.conn.cursor()
        cursor.execute(
            "SELECT subgraph_hash, adjacency_matrix, frequency FROM subgraphs WHERE size = ?",
            (size,),
        )

        results = []
        for row in cursor.fetchall():
            results.append(
                {
                    "hash": row[0],
                    "adjacency_matrix": json.loads(row[1]),
                    "frequency": row[2],
                },
            )

        return results

    def close(self):
        """Close database connection"""
        if self.conn:
            self.conn.close()


class MolecularGraphAnalyser:
    """Analyses molecular graphs and generates subgraphs"""

    def __init__(self, db: SubgraphDatabase):
        self.db = db

    def mol_to_networkx(self, mol: Chem.Mol) -> nx.Graph:
        """Convert RDKit molecule to NetworkX graph"""
        G = nx.Graph()

        # Add nodes (atoms)
        for atom in mol.GetAtoms():
            G.add_node(atom.GetIdx(), symbol=atom.GetSymbol())

        # Add edges (bonds)
        for bond in mol.GetBonds():
            G.add_edge(
                bond.GetBeginAtomIdx(),
                bond.GetEndAtomIdx(),
                bond_type=bond.GetBondType(),
            )

        return G

    def generate_all_subgraphs(
        self,
        G: nx.Graph,
        max_size: int = None,
    ) -> List[nx.Graph]:
        """Generate all connected subgraphs up to max_size, excluding single nodes"""
        if max_size is None:
            max_size = len(G.nodes())

        subgraphs = []
        nodes = list(G.nodes())
        seen_hashes = set()  # Track canonical hashes to avoid duplicates

        # Generate subgraphs of all sizes (starting from 2 nodes to exclude single nodes)
        for size in range(2, min(max_size + 1, len(nodes) + 1)):
            for node_subset in combinations(nodes, size):
                subgraph = G.subgraph(node_subset)
                if nx.is_connected(subgraph) and len(subgraph.edges()) > 0:
                    # Check if we've seen this structure before
                    canonical_hash = self.graph_to_canonical_hash(subgraph)
                    if canonical_hash not in seen_hashes:
                        seen_hashes.add(canonical_hash)
                        subgraphs.append(subgraph.copy())

        return subgraphs

    def graph_to_canonical_hash(self, G: nx.Graph) -> str:
        """Generate canonical hash for a graph using graph isomorphism"""
        # Use NetworkX's weisfeiler_lehman_graph_hash for canonical representation
        # This properly handles graph isomorphism (same structure, different node labels)
        try:
            return nx.weisfeiler_lehman_graph_hash(G)
        except:
            # Fallback to adjacency matrix approach if WL hash fails
            nodes = sorted(G.nodes())
            adj_matrix = []
            for i in nodes:
                row = []
                for j in nodes:
                    if G.has_edge(i, j):
                        row.append(1)
                    else:
                        row.append(0)
                adj_matrix.append(row)
            matrix_str = str(adj_matrix)
            return str(hash(matrix_str))

    def populate_database(self, molecule_input: MoleculeInput):
        """Populate database with subgraphs from a molecule"""
        mol = Chem.MolFromSmiles(molecule_input.smiles)
        G = self.mol_to_networkx(mol)

        logger.info(
            f"Generating subgraphs for {molecule_input.name or molecule_input.smiles}",
        )

        subgraphs = self.generate_all_subgraphs(G)

        for subgraph in subgraphs:
            # Create adjacency matrix for storage
            nodes = sorted(subgraph.nodes())
            adj_matrix = []
            for i in range(len(nodes)):
                row = []
                for j in range(len(nodes)):
                    if subgraph.has_edge(nodes[i], nodes[j]):
                        row.append(1)
                    else:
                        row.append(0)
                adj_matrix.append(row)

            entry = SubgraphEntry(
                subgraph_hash=self.graph_to_canonical_hash(subgraph),
                size=len(subgraph.nodes()),
                adjacency_matrix=adj_matrix,
                parent_molecule=molecule_input.name or molecule_input.smiles,
            )

            self.db.add_subgraph(entry)

        logger.info(f"Added {len(subgraphs)} subgraphs to database")


class GraphEditDistanceCalculator:
    """Calculates graph edit distance using subgraph database"""

    def __init__(self, db: SubgraphDatabase):
        self.db = db

    def calculate_ged_approximation(
        self,
        mol1: MoleculeInput,
        mol2: MoleculeInput,
    ) -> Dict:
        """Calculate approximate graph edit distance using subgraph overlap
        Based on the principle that similar molecules share more subgraphs
        """
        # Convert molecules to graphs
        rdkit_mol1 = Chem.MolFromSmiles(mol1.smiles)
        rdkit_mol2 = Chem.MolFromSmiles(mol2.smiles)

        analyser = MolecularGraphAnalyser(self.db)
        G1 = analyser.mol_to_networkx(rdkit_mol1)
        G2 = analyser.mol_to_networkx(rdkit_mol2)

        # Generate subgraphs for both molecules
        subgraphs1 = analyser.generate_all_subgraphs(G1)
        subgraphs2 = analyser.generate_all_subgraphs(G2)

        # Convert to hash sets for comparison
        hashes1 = {analyser.graph_to_canonical_hash(sg) for sg in subgraphs1}
        hashes2 = {analyser.graph_to_canonical_hash(sg) for sg in subgraphs2}

        # Calculate overlap metrics
        intersection = hashes1.intersection(hashes2)
        union = hashes1.union(hashes2)

        # Tanimoto coefficient (Jaccard similarity)
        tanimoto = len(intersection) / len(union) if union else 0

        # Approximate GED based on subgraph overlap
        # This is a heuristic: fewer shared subgraphs suggest higher edit distance
        max_nodes = max(len(G1.nodes()), len(G2.nodes()))
        min_shared = len(intersection)
        max_possible_shared = min(len(hashes1), len(hashes2))

        # Normalised edit distance approximation
        if max_possible_shared == 0:
            ged_approx = max_nodes
        else:
            similarity_ratio = min_shared / max_possible_shared
            ged_approx = max_nodes * (1 - similarity_ratio)

        return {
            "approximate_ged": round(ged_approx, 3),
            "tanimoto_coefficient": round(tanimoto, 3),
            "shared_subgraphs": len(intersection),
            "total_subgraphs_mol1": len(hashes1),
            "total_subgraphs_mol2": len(hashes2),
            "exact_ged": self._calculate_exact_ged(G1, G2),
        }

    def _calculate_exact_ged(self, G1: nx.Graph, G2: nx.Graph) -> float:
        """Calculate exact graph edit distance using NetworkX"""
        try:
            # This can be computationally expensive for larger graphs
            ged = nx.graph_edit_distance(G1, G2)
            return round(ged, 3)
        except Exception as e:
            logger.warning(f"Could not calculate exact GED: {e}")
            return None


class MolecularVisualiser:
    """Handles molecular visualisation using RDKit"""

    @staticmethod
    def draw_molecule(smiles: str, title: str = "") -> None:
        """Draw a molecule from SMILES string"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"Could not parse SMILES: {smiles}")
            return

        # Generate 2D coordinates
        from rdkit.Chem import rdDepictor

        rdDepictor.Compute2DCoords(mol)

        # Create image
        img = Draw.MolToImage(mol, size=(300, 300))

        # Display using matplotlib
        plt.figure(figsize=(4, 4))
        plt.imshow(img)
        plt.axis("off")
        plt.title(title)
        plt.tight_layout()
        plt.show()

    @staticmethod
    def networkx_to_rdkit_subgraph(
        parent_mol: Chem.Mol,
        subgraph_nodes: List[int],
    ) -> Chem.Mol:
        """Convert NetworkX subgraph back to RDKit mol for visualisation"""
        # Create a new molecule with only the subgraph atoms
        new_mol = Chem.RWMol()

        # Map old atom indices to new ones
        atom_map = {}
        for i, old_idx in enumerate(sorted(subgraph_nodes)):
            atom = parent_mol.GetAtomWithIdx(old_idx)
            new_atom = Chem.Atom(atom.GetSymbol())
            new_idx = new_mol.AddAtom(new_atom)
            atom_map[old_idx] = new_idx

        # Add bonds between atoms in the subgraph
        for bond in parent_mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()

            if begin_idx in subgraph_nodes and end_idx in subgraph_nodes:
                new_mol.AddBond(
                    atom_map[begin_idx],
                    atom_map[end_idx],
                    bond.GetBondType(),
                )

        # Sanitise and return
        try:
            Chem.SanitizeMol(new_mol)
            return new_mol
        except:
            return None

    @staticmethod
    def draw_subgraphs_comparison(
        mol1_smiles: str,
        mol2_smiles: str,
        results: Dict,
        mol1_name: str = "Molecule 1",
        mol2_name: str = "Molecule 2",
    ) -> None:
        """Create comprehensive visualisation with subgraphs ordered by edge count"""
        mol1 = Chem.MolFromSmiles(mol1_smiles)
        mol2 = Chem.MolFromSmiles(mol2_smiles)

        if mol1 is None or mol2 is None:
            logger.error("Could not parse one or both SMILES strings")
            return

        # Generate subgraph data
        db = SubgraphDatabase("molecular_subgraphs.db")
        analyser = MolecularGraphAnalyser(db)

        G1 = analyser.mol_to_networkx(mol1)
        G2 = analyser.mol_to_networkx(mol2)

        subgraphs1 = analyser.generate_all_subgraphs(G1)
        subgraphs2 = analyser.generate_all_subgraphs(G2)

        # Sort subgraphs by edge count (descending), then by node count
        subgraphs1.sort(key=lambda sg: (len(sg.edges()), len(sg.nodes())), reverse=True)
        subgraphs2.sort(key=lambda sg: (len(sg.edges()), len(sg.nodes())), reverse=True)

        # Group subgraphs by edge count for display
        sg1_by_edges = {}
        sg2_by_edges = {}

        for sg in subgraphs1:
            edge_count = len(sg.edges())
            if edge_count not in sg1_by_edges:
                sg1_by_edges[edge_count] = []
            sg1_by_edges[edge_count].append(sg)

        for sg in subgraphs2:
            edge_count = len(sg.edges())
            if edge_count not in sg2_by_edges:
                sg2_by_edges[edge_count] = []
            sg2_by_edges[edge_count].append(sg)

        # Find shared subgraphs using canonical hashes
        hashes1 = {analyser.graph_to_canonical_hash(sg): sg for sg in subgraphs1}
        hashes2 = {analyser.graph_to_canonical_hash(sg): sg for sg in subgraphs2}
        shared_hashes = set(hashes1.keys()).intersection(set(hashes2.keys()))

        # Create main comparison figure
        fig = plt.figure(figsize=(20, 14))

        # Main molecules
        from rdkit.Chem import rdDepictor

        rdDepictor.Compute2DCoords(mol1)
        rdDepictor.Compute2DCoords(mol2)

        # Top row: main molecules
        ax_mol1 = plt.subplot2grid((5, 8), (0, 1), colspan=2)
        ax_mol2 = plt.subplot2grid((5, 8), (0, 5), colspan=2)

        img1 = Draw.MolToImage(mol1, size=(300, 300))
        img2 = Draw.MolToImage(mol2, size=(300, 300))

        ax_mol1.imshow(img1)
        ax_mol1.set_title(f"{mol1_name}\n{mol1_smiles}", fontsize=12, fontweight="bold")
        ax_mol1.axis("off")

        ax_mol2.imshow(img2)
        ax_mol2.set_title(f"{mol2_name}\n{mol2_smiles}", fontsize=12, fontweight="bold")
        ax_mol2.axis("off")

        # Results summary in the middle
        ax_results = plt.subplot2grid((5, 8), (0, 3), colspan=2)

        # Count unique subgraphs by edge count
        edges1_counts = {k: len(v) for k, v in sg1_by_edges.items()}
        edges2_counts = {k: len(v) for k, v in sg2_by_edges.items()}

        results_text = f"""GED Analysis Results
        
Approximate GED: {results["approximate_ged"]}
Exact GED: {results.get("exact_ged", "N/A")}
Tanimoto Coefficient: {results["tanimoto_coefficient"]}

Shared Subgraphs: {results["shared_subgraphs"]}
Unique subgraphs in {mol1_name}: {len(subgraphs1)}
Unique subgraphs in {mol2_name}: {len(subgraphs2)}

Edge distribution:
{mol1_name}: {dict(sorted(edges1_counts.items(), reverse=True))}
{mol2_name}: {dict(sorted(edges2_counts.items(), reverse=True))}"""

        ax_results.text(
            0.5,
            0.5,
            results_text,
            ha="center",
            va="center",
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.5", facecolor="lightblue", alpha=0.7),
        )
        ax_results.axis("off")

        # Subgraphs visualisation ordered by edge count (descending)
        all_edge_counts = sorted(
            set(list(sg1_by_edges.keys()) + list(sg2_by_edges.keys())),
            reverse=True,
        )

        row_offset = 1
        for edge_count in all_edge_counts:
            if row_offset >= 5:  # Prevent overflow
                break

            # Edge count label
            ax_size = plt.subplot2grid((5, 8), (row_offset, 0))
            ax_size.text(
                0.5,
                0.5,
                f"{edge_count}\nedges",
                ha="center",
                va="center",
                fontsize=11,
                fontweight="bold",
            )
            ax_size.axis("off")

            # Molecule 1 subgraphs for this edge count
            mol1_sgs = sg1_by_edges.get(edge_count, [])
            col_start = 1
            for i, sg in enumerate(mol1_sgs[:3]):  # Limit to 3 per row
                if col_start + i < 4:
                    ax = plt.subplot2grid((5, 8), (row_offset, col_start + i))

                    # Check if this subgraph is shared
                    sg_hash = analyser.graph_to_canonical_hash(sg)
                    is_shared = sg_hash in shared_hashes

                    # Convert subgraph to RDKit mol for drawing
                    sg_mol = MolecularVisualiser.networkx_to_rdkit_subgraph(
                        mol1,
                        list(sg.nodes()),
                    )

                    if sg_mol is not None:
                        rdDepictor.Compute2DCoords(sg_mol)
                        sg_img = Draw.MolToImage(sg_mol, size=(120, 120))
                        ax.imshow(sg_img)

                        # Colour border based on sharing
                        border_colour = "red" if is_shared else "lightgray"
                        border_width = 3 if is_shared else 1

                        for spine in ax.spines.values():
                            spine.set_edgecolor(border_colour)
                            spine.set_linewidth(border_width)

                    # Show structure info
                    nodes_sorted = sorted(list(sg.nodes()))
                    edges_list = [(min(e), max(e)) for e in sg.edges()]
                    edges_sorted = sorted(edges_list)
                    ax.set_title(f"N:{nodes_sorted}\nE:{edges_sorted}", fontsize=7)
                    ax.axis("off")

            # Molecule 2 subgraphs for this edge count
            mol2_sgs = sg2_by_edges.get(edge_count, [])
            col_start = 5
            for i, sg in enumerate(mol2_sgs[:3]):  # Limit to 3 per row
                if col_start + i < 8:
                    ax = plt.subplot2grid((5, 8), (row_offset, col_start + i))

                    # Check if this subgraph is shared
                    sg_hash = analyser.graph_to_canonical_hash(sg)
                    is_shared = sg_hash in shared_hashes

                    # Convert subgraph to RDKit mol for drawing
                    sg_mol = MolecularVisualiser.networkx_to_rdkit_subgraph(
                        mol2,
                        list(sg.nodes()),
                    )

                    if sg_mol is not None:
                        rdDepictor.Compute2DCoords(sg_mol)
                        sg_img = Draw.MolToImage(sg_mol, size=(120, 120))
                        ax.imshow(sg_img)

                        # Colour border based on sharing
                        border_colour = "red" if is_shared else "lightgray"
                        border_width = 3 if is_shared else 1

                        for spine in ax.spines.values():
                            spine.set_edgecolor(border_colour)
                            spine.set_linewidth(border_width)

                    # Show structure info
                    nodes_sorted = sorted(list(sg.nodes()))
                    edges_list = [(min(e), max(e)) for e in sg.edges()]
                    edges_sorted = sorted(edges_list)
                    ax.set_title(f"N:{nodes_sorted}\nE:{edges_sorted}", fontsize=7)
                    ax.axis("off")

            row_offset += 1

        # Legend
        legend_ax = plt.subplot2grid((5, 8), (4, 3), colspan=2)
        legend_text = """Red border = Shared subgraph (same structure)
Gray border = Unique subgraph
N: = Node indices, E: = Edge pairs
Ordered by edge count (descending)"""
        legend_ax.text(
            0.5,
            0.5,
            legend_text,
            ha="center",
            va="center",
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow"),
        )
        legend_ax.axis("off")

        plt.suptitle(
            f"Molecular Subgraph Comparison (by Edge Count): {mol1_name} vs {mol2_name}",
            fontsize=16,
            fontweight="bold",
        )
        plt.tight_layout()
        plt.show()

        # Print detailed subgraph analysis
        print("\nDetailed Subgraph Analysis:")
        print(f"{'=' * 50}")
        print(f"{mol1_name} ({mol1_smiles}):")
        for edge_count in sorted(sg1_by_edges.keys(), reverse=True):
            sgs = sg1_by_edges[edge_count]
            print(f"  {edge_count} edges: {len(sgs)} unique subgraph(s)")
            for sg in sgs[:3]:  # Show first 3
                nodes = sorted(list(sg.nodes()))
                edges = sorted([(min(e), max(e)) for e in sg.edges()])
                is_shared = analyser.graph_to_canonical_hash(sg) in shared_hashes
                shared_mark = " [SHARED]" if is_shared else ""
                print(f"    Nodes: {nodes}, Edges: {edges}{shared_mark}")

        print(f"\n{mol2_name} ({mol2_smiles}):")
        for edge_count in sorted(sg2_by_edges.keys(), reverse=True):
            sgs = sg2_by_edges[edge_count]
            print(f"  {edge_count} edges: {len(sgs)} unique subgraph(s)")
            for sg in sgs[:3]:  # Show first 3
                nodes = sorted(list(sg.nodes()))
                edges = sorted([(min(e), max(e)) for e in sg.edges()])
                is_shared = analyser.graph_to_canonical_hash(sg) in shared_hashes
                shared_mark = " [SHARED]" if is_shared else ""
                print(f"    Nodes: {nodes}, Edges: {edges}{shared_mark}")

        db.close()

    @staticmethod
    def compare_molecules(
        mol1_smiles: str,
        mol2_smiles: str,
        results: Dict,
        mol1_name: str = "Molecule 1",
        mol2_name: str = "Molecule 2",
    ) -> None:
        """Create side-by-side comparison of two molecules with results"""
        # Use the enhanced subgraph comparison
        MolecularVisualiser.draw_subgraphs_comparison(
            mol1_smiles,
            mol2_smiles,
            results,
            mol1_name,
            mol2_name,
        )


def main():
    """Main execution function with example usage"""
    # Initialise database and components
    db = SubgraphDatabase("molecular_subgraphs.db")
    analyser = MolecularGraphAnalyser(db)
    calculator = GraphEditDistanceCalculator(db)
    visualiser = MolecularVisualiser()

    try:
        # Example molecules (carbon chains and simple structures)
        examples = [
            MoleculeInput(smiles="CCCC", name="Butane"),
            MoleculeInput(smiles="CCC(C)C", name="Isopentane"),
            MoleculeInput(smiles="C1CCC1", name="Cyclobutane"),
            MoleculeInput(smiles="CC(C)(C)C", name="Neopentane"),
            MoleculeInput(smiles="CCCCCC", name="Hexane"),
        ]

        # Populate database with subgraphs
        logger.info("Populating subgraph database...")
        for mol in examples:
            analyser.populate_database(mol)

        # Calculate graph edit distances between pairs
        logger.info("Calculating graph edit distances...")

        # Example comparison: Butane vs Isopentane
        mol1 = examples[3]  # Butane
        mol2 = examples[4]  # Isopentane

        results = calculator.calculate_ged_approximation(mol1, mol2)

        print("\nGraph Edit Distance Analysis:")
        print(f"Molecule 1: {mol1.name} ({mol1.smiles})")
        print(f"Molecule 2: {mol2.name} ({mol2.smiles})")
        print(f"Approximate GED: {results['approximate_ged']}")
        print(f"Exact GED: {results.get('exact_ged', 'Not calculated')}")
        print(f"Tanimoto Coefficient: {results['tanimoto_coefficient']}")
        print(f"Shared Subgraphs: {results['shared_subgraphs']}")

        # Visualise comparison
        visualiser.compare_molecules(
            mol1.smiles,
            mol2.smiles,
            results,
            mol1.name,
            mol2.name,
        )

        # Additional comparisons
        print("\n" + "=" * 50)
        print("Additional Comparisons:")
        print("=" * 50)

        for i in range(len(examples)):
            for j in range(i + 1, len(examples)):
                mol_a, mol_b = examples[i], examples[j]
                results = calculator.calculate_ged_approximation(mol_a, mol_b)
                print(
                    f"{mol_a.name} vs {mol_b.name}: GED ≈ {results['approximate_ged']}, "
                    f"Tanimoto = {results['tanimoto_coefficient']}",
                )

    except Exception as e:
        logger.error(f"Error in main execution: {e}")
        raise
    finally:
        db.close()


if __name__ == "__main__":
    main()
