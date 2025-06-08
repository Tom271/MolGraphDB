"""FastAPI backend for MolGraphDB - Molecular Similarity Calculator.

This backend provides API endpoints for calculating molecular similarity using graph
edit distance.
"""

import base64
import io
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field, validator

# Try to import RDKit for molecule visualization
try:
    from rdkit import Chem
    from rdkit.Chem import Draw, rdDepictor

    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logger = logging.getLogger(__name__)
    logger.warning("RDKit not available - molecule images will not be generated")

# Add the src directory to Python path to import mcs module
src_path = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(src_path))

# Import the calculation functions from mcs.py
try:
    from mcs import (
        GraphEditDistanceCalculator,
        MolecularGraphAnalyser,
        MoleculeInput,
        SubgraphDatabase,
    )
    MCS_AVAILABLE = True
except ImportError as e:
    logging.exception(f"Failed to import mcs module: {e}")
    MCS_AVAILABLE = False

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Create FastAPI app
app = FastAPI(
    title="MolGraphDB API",
    description="Molecular similarity calculator using graph edit distance",
    version="1.0.0",
)

# Enable CORS for the React frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000"],  # React dev server
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# Request/Response models
class MolecularSimilarityRequest(BaseModel):
    smiles1: str = Field(..., description="First molecule SMILES string")
    smiles2: str = Field(..., description="Second molecule SMILES string")

    @validator("smiles1", "smiles2")
    def validate_smiles(cls, v):
        if not v or not v.strip():
            raise ValueError("SMILES string cannot be empty")
        return v.strip()


class SubgraphVisualization(BaseModel):
    id: str = Field(..., description="Unique identifier for the subgraph")
    nodes: List[int] = Field(..., description="Node indices in the subgraph")
    edges: List[List[int]] = Field(..., description="Edge connections")
    edge_count: int = Field(..., description="Number of edges in subgraph")
    node_count: int = Field(..., description="Number of nodes in subgraph")
    is_shared: bool = Field(
        ...,
        description="Whether this subgraph is shared between molecules",
    )
    smiles: Optional[str] = Field(None, description="SMILES representation of subgraph")
    image: Optional[str] = Field(None, description="Base64 encoded image of subgraph")


class MoleculeVisualization(BaseModel):
    smiles: str = Field(..., description="SMILES string")
    image: Optional[str] = Field(None, description="Base64 encoded molecule image")
    subgraphs: List[SubgraphVisualization] = Field(
        ...,
        description="List of subgraphs organized by edge count",
    )


class MolecularSimilarityResponse(BaseModel):
    similarity_score: float = Field(..., description="Similarity score (0-1)")
    graph_edit_distance: float = Field(..., description="Graph edit distance")
    subgraphs1_count: int = Field(..., description="Number of subgraphs in molecule 1")
    subgraphs2_count: int = Field(..., description="Number of subgraphs in molecule 2")
    molecule1_info: Dict[str, Any] = Field(
        ...,
        description="Information about molecule 1",
    )
    molecule2_info: Dict[str, Any] = Field(
        ...,
        description="Information about molecule 2",
    )
    molecule1_viz: MoleculeVisualization = Field(
        ...,
        description="Visualization data for molecule 1",
    )
    molecule2_viz: MoleculeVisualization = Field(
        ...,
        description="Visualization data for molecule 2",
    )
    shared_subgraphs: List[str] = Field(
        ...,
        description="IDs of subgraphs shared between both molecules",
    )


# Global database instance
db = None


def generate_molecule_image(smiles: str, size: tuple = (300, 300)) -> Optional[str]:
    """Generate base64 encoded image of molecule from SMILES"""
    if not RDKIT_AVAILABLE:
        return None

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # Generate 2D coordinates
        rdDepictor.Compute2DCoords(mol)

        # Create image
        img = Draw.MolToImage(mol, size=size)

        # Convert to base64
        buffered = io.BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode()

        return f"data:image/png;base64,{img_str}"
    except Exception as e:
        logger.warning(f"Failed to generate molecule image: {e}")
        return None


def generate_subgraph_image(
    parent_mol,
    subgraph_nodes: List[int],
    size: tuple = (150, 150),
) -> Optional[str]:
    """Generate base64 encoded image of subgraph"""
    if not RDKIT_AVAILABLE:
        return None

    try:
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

        # Sanitize and generate image
        try:
            Chem.SanitizeMol(new_mol)
            rdDepictor.Compute2DCoords(new_mol)
            img = Draw.MolToImage(new_mol, size=size)

            # Convert to base64
            buffered = io.BytesIO()
            img.save(buffered, format="PNG")
            img_str = base64.b64encode(buffered.getvalue()).decode()

            return f"data:image/png;base64,{img_str}"
        except:
            return None

    except Exception as e:
        logger.warning(f"Failed to generate subgraph image: {e}")
        return None


def get_subgraph_smiles(parent_mol, subgraph_nodes: List[int]) -> Optional[str]:
    """Generate SMILES string for subgraph"""
    if not RDKIT_AVAILABLE:
        return None

    try:
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

        # Sanitize and generate SMILES
        try:
            Chem.SanitizeMol(new_mol)
            return Chem.MolToSmiles(new_mol)
        except:
            return None

    except Exception as e:
        logger.warning(f"Failed to generate subgraph SMILES: {e}")
        return None



def get_molecule_info(smiles: str) -> dict[str, Any]:
    """Get comprehensive molecular information from SMILES string"""
    if not RDKIT_AVAILABLE:
        return {
            "smiles": smiles,
            "molecular_formula": "Unknown",
            "molecular_weight": 0.0,
            "num_atoms": 0,
            "num_bonds": 0,
        }
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {
                "smiles": smiles,
                "molecular_formula": "Invalid",
                "molecular_weight": 0.0,
                "num_atoms": 0,
                "num_bonds": 0,
            }
        
        # Import additional RDKit modules for descriptors
        from rdkit.Chem import rdMolDescriptors
        
        # Count atoms by type
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "C")
        
        return {
            "smiles": smiles,
            "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
            "molecular_weight": round(rdMolDescriptors.CalcExactMolWt(mol), 2),
            "num_atoms": mol.GetNumAtoms(),
            "num_bonds": mol.GetNumBonds(),
            "carbon_atoms": carbon_count,
        }
    
    except Exception as e:
        logger.warning(f"Failed to get molecule info for {smiles}: {e}")
        return {
            "smiles": smiles,
            "molecular_formula": "Error",
            "molecular_weight": 0.0,
            "num_atoms": 0,
            "num_bonds": 0,
        }

def initialize_database():
    """Initialize the subgraph database"""
    global db
    try:
        db_path = Path(__file__).parent.parent / "molecular_subgraphs.db"
        if MCS_AVAILABLE:
            db = SubgraphDatabase(str(db_path))
            logger.info(f"Database initialized at {db_path}")
        else:
            logger.warning("MCS module not available, database will not be initialized")
            db = None
    except Exception as e:
        logger.error(f"Failed to initialize database: {e}")
        db = None


def create_visualization_data(
    smiles: str,
    subgraphs_data: List,
    shared_subgraph_ids: set,
) -> MoleculeVisualization:
    """Create visualization data for a molecule including subgraphs organized by edge count"""
    if not RDKIT_AVAILABLE:
        return MoleculeVisualization(smiles=smiles, image=None, subgraphs=[])

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return MoleculeVisualization(smiles=smiles, image=None, subgraphs=[])

    # Generate molecule image
    mol_image = generate_molecule_image(smiles)

    # Create subgraph visualizations
    subgraph_viz = []
    for sg_data in subgraphs_data:
        # Generate subgraph image and SMILES
        sg_image = generate_subgraph_image(mol, sg_data["nodes"])
        sg_smiles = get_subgraph_smiles(mol, sg_data["nodes"])

        # Create edges list from node connections
        edges = []
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if begin_idx in sg_data["nodes"] and end_idx in sg_data["nodes"]:
                edges.append([begin_idx, end_idx])

        subgraph_viz.append(
            SubgraphVisualization(
                id=sg_data["id"],
                nodes=sg_data["nodes"],
                edges=edges,
                edge_count=len(edges),
                node_count=len(sg_data["nodes"]),
                is_shared=sg_data["id"] in shared_subgraph_ids,
                smiles=sg_smiles,
                image=sg_image,
            ),
        )

    # Sort subgraphs by edge count (descending) then by node count
    subgraph_viz.sort(key=lambda x: (x.edge_count, x.node_count), reverse=True)

    return MoleculeVisualization(smiles=smiles, image=mol_image, subgraphs=subgraph_viz)


def calculate_molecular_similarity_mcs(smiles1: str, smiles2: str) -> dict[str, Any]:
    """Calculate molecular similarity using the mcs.py module
    This function calls the actual implementation from mcs.py
    """
    # Get molecule information
    mol1_info = get_molecule_info(smiles1)
    mol2_info = get_molecule_info(smiles2)
    
    # Check if we have everything needed for full MCS calculation
    if not MCS_AVAILABLE or db is None:
        logger.info("MCS module or database not available, using placeholder calculation")
        return calculate_molecular_similarity_placeholder(smiles1, smiles2)
    
    try:
        # Create molecule inputs
        mol1 = MoleculeInput(smiles=smiles1, name="Molecule 1")
        mol2 = MoleculeInput(smiles=smiles2, name="Molecule 2")

        # Initialize components
        analyser = MolecularGraphAnalyser(db)
        calculator = GraphEditDistanceCalculator(db)

        # Calculate similarity
        results = calculator.calculate_ged_approximation(mol1, mol2)

        # Generate subgraphs for visualization if RDKit is available
        if RDKIT_AVAILABLE:
            mol1_rdkit = Chem.MolFromSmiles(smiles1)
            mol2_rdkit = Chem.MolFromSmiles(smiles2)

            if mol1_rdkit and mol2_rdkit:
                # Get subgraphs from the analyser
                g1 = analyser.mol_to_networkx(mol1_rdkit)
                g2 = analyser.mol_to_networkx(mol2_rdkit)

                subgraphs1 = analyser.generate_all_subgraphs(g1)
                subgraphs2 = analyser.generate_all_subgraphs(g2)

                # Create hash mapping for shared subgraphs
                hashes1 = {
                    analyser.graph_to_canonical_hash(sg): sg for sg in subgraphs1
                }
                hashes2 = {
                    analyser.graph_to_canonical_hash(sg): sg for sg in subgraphs2
                }
                shared_hashes = set(hashes1.keys()).intersection(
                    set(hashes2.keys()),
                )

                # Convert to visualization format
                subgraphs1_data = [
                    {
                        "id": analyser.graph_to_canonical_hash(sg),
                        "nodes": list(sg.nodes()),
                    }
                    for sg in subgraphs1
                ]

                subgraphs2_data = [
                    {
                        "id": analyser.graph_to_canonical_hash(sg),
                        "nodes": list(sg.nodes()),
                    }
                    for sg in subgraphs2
                ]

                # Create visualization data
                mol1_viz = create_visualization_data(
                    smiles1,
                    subgraphs1_data,
                    shared_hashes,
                )
                mol2_viz = create_visualization_data(
                    smiles2,
                    subgraphs2_data,
                    shared_hashes,
                )
            else:
                mol1_viz = MoleculeVisualization(
                    smiles=smiles1,
                    image=generate_molecule_image(smiles1),
                    subgraphs=[],
                )
                mol2_viz = MoleculeVisualization(
                    smiles=smiles2,
                    image=generate_molecule_image(smiles2),
                    subgraphs=[],
                )
                shared_hashes = set()
        else:
            # No RDKit available, create empty visualization
            mol1_viz = MoleculeVisualization(
                smiles=smiles1,
                image=None,
                subgraphs=[],
            )
            mol2_viz = MoleculeVisualization(
                smiles=smiles2,
                image=None,
                subgraphs=[],
            )
            shared_hashes = set()

        return {
            "similarity_score": 1.0
            - (results["approximate_ged"] / 10.0),  # Convert to similarity
            "graph_edit_distance": results["approximate_ged"],
            "subgraphs1_count": results["total_subgraphs_mol1"],
            "subgraphs2_count": results["total_subgraphs_mol2"],
            "molecule1_info": mol1_info,
            "molecule2_info": mol2_info,
            "molecule1_viz": mol1_viz,
            "molecule2_viz": mol2_viz,
            "shared_subgraphs": list(shared_hashes),
        }

    except Exception as e:
        logger.error(f"Error in mcs calculation: {e}")
        # Fallback to placeholder
        return calculate_molecular_similarity_placeholder(smiles1, smiles2)


def calculate_molecular_similarity_placeholder(
    smiles1: str,
    smiles2: str,
) -> dict[str, Any]:
    """Placeholder function for molecular similarity calculation
    This will be replaced with actual calculations from mcs.py
    """
    logger.info(f"Using placeholder calculation for {smiles1} vs {smiles2}")

    # Simple placeholder logic based on SMILES length difference
    len_diff = abs(len(smiles1) - len(smiles2))
    max_len = max(len(smiles1), len(smiles2))

    # Mock similarity score (higher for similar lengths)
    similarity = max(0.0, 1.0 - (len_diff / max_len)) if max_len > 0 else 1.0

    # Mock graph edit distance
    ged = len_diff * 1.5

    # Get actual molecule info using our function
    mol1_info = get_molecule_info(smiles1)
    mol2_info = get_molecule_info(smiles2)

    # Create placeholder visualization data
    mol1_viz = MoleculeVisualization(
        smiles=smiles1,
        image=generate_molecule_image(smiles1) if RDKIT_AVAILABLE else None,
        subgraphs=[],
    )
    mol2_viz = MoleculeVisualization(
        smiles=smiles2,
        image=generate_molecule_image(smiles2) if RDKIT_AVAILABLE else None,
        subgraphs=[],
    )

    return {
        "similarity_score": similarity,
        "graph_edit_distance": ged,
        "subgraphs1_count": max(1, len(smiles1) // 3),  # Mock subgraph count
        "subgraphs2_count": max(1, len(smiles2) // 3),  # Mock subgraph count
        "molecule1_info": mol1_info,
        "molecule2_info": mol2_info,
        "molecule1_viz": mol1_viz,
        "molecule2_viz": mol2_viz,
        "shared_subgraphs": [],
    }


@app.on_event("startup")
async def startup_event():
    """Initialize database on startup"""
    initialize_database()


@app.get("/")
async def root():
    """Root endpoint"""
    return {
        "message": "MolGraphDB API - Molecular Similarity Calculator",
        "version": "1.0.0",
        "docs": "/docs",
    }


@app.get("/health")
async def health_check():
    """Health check endpoint"""
    return {"status": "healthy", "database_connected": db is not None}


@app.post("/api/calculate_similarity", response_model=MolecularSimilarityResponse)
async def calculate_similarity(request: MolecularSimilarityRequest):
    """Calculate molecular similarity between two SMILES strings"""
    try:
        logger.info(f"Calculating similarity: {request.smiles1} vs {request.smiles2}")

        # Call the calculation function (this will use mcs.py when available)
        results = calculate_molecular_similarity_mcs(request.smiles1, request.smiles2)

        return MolecularSimilarityResponse(**results)

    except ValueError as e:
        logger.error(f"Validation error: {e}")
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Calculation error: {e}")
        raise HTTPException(status_code=500, detail=f"Calculation failed: {e!s}")


if __name__ == "__main__":
    import uvicorn

    # Run the server
    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=True, log_level="info")
