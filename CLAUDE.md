# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

### Python Development (UV package manager)
- **Run tests**: `uv run pytest` or `uv run pytest tests/`
- **Run specific test**: `uv run pytest tests/test_hello.py`
- **Install dependencies**: `uv sync`
- **Lint code**: `uv run ruff check src/ backend/`
- **Format code**: `uv run ruff format src/ backend/`
- **Start backend API**: `cd backend && uv run python main.py`
- **Run main MCS script**: `uv run python src/mcs.py`

### Frontend Development (Node.js/npm)
- **Install dependencies**: `cd frontend && npm install`
- **Start development server**: `cd frontend && npm start`
- **Build for production**: `cd frontend && npm run build`
- **Run tests**: `cd frontend && npm test`

### Julia Development
- **Run Julia code**: `julia juliaMolGraph/src.jl`
- **Julia REPL with project**: `julia --project=juliaMolGraph`

### Database
- **SQLite database file**: `molecular_subgraphs.db` (contains molecular subgraph data)

## Project Architecture

### Core Components

**Molecular Graph Analysis (`src/mcs.py`)**
- Main molecular graph edit distance calculator
- Uses NetworkX for graph operations and RDKit for molecular parsing
- Implements subgraph-based similarity calculation inspired by NextMove SmallWorld
- Database-backed subgraph storage using SQLite
- Comprehensive visualization system with matplotlib

**Backend API (`backend/main.py`)**
- FastAPI web service for molecular similarity calculations
- RESTful API endpoints for frontend integration
- CORS enabled for React frontend at localhost:3000
- Placeholder calculations with fallback to actual mcs.py implementation
- Pydantic models for request/response validation

**Frontend (`frontend/src/App.js`)**
- React component for molecular similarity calculator UI
- Clean, responsive design with custom CSS
- SMILES input validation and error handling
- Real-time API integration with loading states and error handling
- Visual display of similarity scores and molecular information

**Parallel Julia Implementation (`juliaMolGraph/`)**
- Julia version using MolecularGraph.jl, Graphs.jl
- Same core functionality as Python implementation
- Uses DataFrames for data manipulation and SQLite for storage

### Key Design Patterns

**Molecular Input Validation**
- Pydantic models with custom validators
- Strict constraints: carbon-only molecules, single bonds, max 20 atoms
- SMILES string parsing and sanitization

**Subgraph Generation Strategy**
- Exhaustive connected subgraph enumeration (size 2+ nodes)
- Canonical hash generation using Weisfeiler-Lehman algorithm
- Database deduplication with frequency tracking

**Graph Edit Distance Calculation**
- Exact GED using NetworkX with atom/bond type matching
- Subgraph-based approximation for efficiency
- Tanimoto coefficient calculation for similarity scoring

### Database Schema
- **subgraphs table**: stores canonical subgraph representations
- **Fields**: hash, size, adjacency_matrix, parent_molecules, frequency
- **Indexes**: on hash and size for query optimization

### Visualization System
- RDKit molecule image generation (PNG format)
- NetworkX subgraph visualization with shared structure highlighting
- Comprehensive comparison plots showing edge count distribution

## Development Notes

**Dependencies**
- Python: RDKit, NetworkX, FastAPI, Pydantic, matplotlib, SQLite3, uvicorn
- Julia: MolecularGraph.jl, Graphs.jl, DataFrames.jl, SQLite.jl
- Frontend: React, react-scripts (Create React App)

**Performance Considerations**
- Subgraph generation is computationally expensive (exponential in molecule size)
- Database caching reduces redundant calculations
- Frontend limits subgraph display to first 50 for performance
- Graph edit distance calculation has 30-second timeout

**Testing**
- Use pytest for Python tests
- Test coverage includes hello.py functionality
- Molecular constraint validation should be tested with edge cases