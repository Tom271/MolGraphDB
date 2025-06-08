#!/usr/bin/env julia
"""
Molecular Graph Edit Distance Calculator
Inspired by NextMove SmallWorld approach for subgraph-based similarity

This project calculates graph edit distance between carbon-only molecules
using a database of molecular subgraphs for efficient computation.
"""

using SQLite
using DataFrames
using Graphs
using MetaGraphs
using JSON3
using SHA
using Logging
using Combinatorics
using MolecularGraph
using GraphsMatching
using Plots
using StatsBase
using LinearAlgebra

# Configure logging
logger = ConsoleLogger()
global_logger(logger)

# Custom structures for molecular data
struct MoleculeInput
    smiles::String
    name::Union{String,Nothing}

    function MoleculeInput(smiles::String, name::Union{String,Nothing}=nothing)
        # Validate SMILES - basic validation for now
        mol = try
            smilestomol(smiles)
        catch e
            throw(ArgumentError("Invalid SMILES string: $smiles"))
        end

        # Check atom count
        if nv(mol) > 15
            throw(ArgumentError("Molecule has $(nv(mol)) atoms, maximum is 15"))
        end

        # Check only carbon atoms (simplified for this example)
        # In practice, you'd iterate through atoms and check element types

        new(smiles, name)
    end
end

struct SubgraphEntry
    subgraph_hash::String
    size::Int
    adjacency_matrix::Matrix{Int}
    parent_molecule::String
    frequency::Int
end

# Database management
mutable struct SubgraphDatabase
    db_path::String
    conn::SQLite.DB

    function SubgraphDatabase(db_path::String="subgraphs.db")
        conn = SQLite.DB(db_path)
        db = new(db_path, conn)
        init_database!(db)
        return db
    end
end

function init_database!(db::SubgraphDatabase)
    """Initialise the SQLite database"""
    SQLite.execute(
        db.conn,
        """
    CREATE TABLE IF NOT EXISTS subgraphs (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        subgraph_hash TEXT UNIQUE,
        size INTEGER,
        adjacency_matrix TEXT,
        parent_molecules TEXT,
        frequency INTEGER DEFAULT 1
    )
"""
    )

    SQLite.execute(db.conn, "CREATE INDEX IF NOT EXISTS idx_hash ON subgraphs(subgraph_hash)")
    SQLite.execute(db.conn, "CREATE INDEX IF NOT EXISTS idx_size ON subgraphs(size)")

    @info "Database initialised"
end

function add_subgraph!(db::SubgraphDatabase, entry::SubgraphEntry)
    """Add or update a subgraph entry"""
    # Check if subgraph already exists
    result = SQLite.execute(db.conn,
        "SELECT id, parent_molecules, frequency FROM subgraphs WHERE subgraph_hash = ?",
        [entry.subgraph_hash])

    if !isempty(result)
        # Update existing entry
        row = first(result)
        existing_id, existing_parents, existing_freq = row.id, row.parent_molecules, row.frequency
        parents_list = JSON3.read(existing_parents, Vector{String})

        if !(entry.parent_molecule in parents_list)
            push!(parents_list, entry.parent_molecule)
        end

        SQLite.execute(db.conn, """
            UPDATE subgraphs 
            SET parent_molecules = ?, frequency = ?
            WHERE id = ?
        """, [JSON3.write(parents_list), existing_freq + 1, existing_id])
    else
        # Insert new entry
        SQLite.execute(db.conn, """
            INSERT INTO subgraphs 
            (subgraph_hash, size, adjacency_matrix, parent_molecules, frequency)
            VALUES (?, ?, ?, ?, ?)
        """, [
                entry.subgraph_hash,
                entry.size,
                JSON3.write(entry.adjacency_matrix),
                JSON3.write([entry.parent_molecule]),
                entry.frequency
            ])
    end
end

function get_subgraphs_by_size(db::SubgraphDatabase, size::Int)
    """Retrieve all subgraphs of a given size"""
    result = SQLite.execute(db.conn,
        "SELECT subgraph_hash, adjacency_matrix, frequency FROM subgraphs WHERE size = ?",
        [size])

    return [(
        hash=row.subgraph_hash,
        adjacency_matrix=JSON3.read(row.adjacency_matrix, Matrix{Int}),
        frequency=row.frequency
    ) for row in result]
end

function Base.close(db::SubgraphDatabase)
    """Close database connection"""
    SQLite.close(db.conn)
end

# Molecular graph analysis
struct MolecularGraphAnalyser
    db::SubgraphDatabase
end

function mol_to_graph(mol_input::MoleculeInput)
    """Convert SMILES to Graphs.jl graph"""
    # Using MolecularGraph.jl for SMILES parsing
    mol = smilestomol(mol_input.smiles)

    # Convert to simple graph
    n_atoms = nv(mol)
    g = SimpleGraph(n_atoms)

    # Add edges based on molecular bonds
    for edge in edges(mol)
        add_edge!(g, src(edge), dst(edge))
    end

    return g
end

function generate_all_subgraphs(g::SimpleGraph, max_size::Union{Int,Nothing}=nothing)
    """Generate all connected subgraphs up to max_size, excluding single nodes"""
    n_nodes = nv(g)
    max_size = isnothing(max_size) ? n_nodes : max_size

    subgraphs = SimpleGraph[]
    nodes = collect(1:n_nodes)
    seen_hashes = Set{String}()

    # Generate subgraphs of all sizes (starting from 2 nodes)
    for size in 2:min(max_size, n_nodes)
        for node_subset in combinations(nodes, size)
            # Create induced subgraph
            subgraph_nodes = collect(node_subset)
            sg = induced_subgraph(g, subgraph_nodes)[1]

            # Check if connected and has edges
            if is_connected(sg) && ne(sg) > 0
                canonical_hash = graph_to_canonical_hash(sg)
                if !(canonical_hash in seen_hashes)
                    push!(seen_hashes, canonical_hash)
                    push!(subgraphs, sg)
                end
            end
        end
    end

    return subgraphs
end

function graph_to_canonical_hash(g::SimpleGraph)
    """Generate canonical hash for a graph"""
    # Use adjacency matrix as a simple canonical representation
    # In practice, you might want a more sophisticated graph canonicalization
    adj_matrix = adjacency_matrix(g)

    # Sort nodes by degree for some canonicalization
    degrees = [degree(g, v) for v in vertices(g)]
    perm = sortperm(degrees, rev=true)

    # Reorder adjacency matrix
    canonical_adj = adj_matrix[perm, perm]

    # Generate hash
    return bytes2hex(sha256(string(canonical_adj)))
end

function populate_database!(analyser::MolecularGraphAnalyser, molecule_input::MoleculeInput)
    """Populate database with subgraphs from a molecule"""
    g = mol_to_graph(molecule_input)

    @info "Generating subgraphs for $(something(molecule_input.name, molecule_input.smiles))"

    subgraphs = generate_all_subgraphs(g)

    for sg in subgraphs
        # Create adjacency matrix
        adj_matrix = Matrix{Int}(adjacency_matrix(sg))

        entry = SubgraphEntry(
            graph_to_canonical_hash(sg),
            nv(sg),
            adj_matrix,
            something(molecule_input.name, molecule_input.smiles),
            1
        )

        add_subgraph!(analyser.db, entry)
    end

    @info "Added $(length(subgraphs)) subgraphs to database"
end

# Graph edit distance calculation
struct GraphEditDistanceCalculator
    db::SubgraphDatabase
end

function calculate_ged_approximation(calc::GraphEditDistanceCalculator,
    mol1::MoleculeInput,
    mol2::MoleculeInput)
    """Calculate approximate graph edit distance using subgraph overlap"""

    # Convert molecules to graphs
    g1 = mol_to_graph(mol1)
    g2 = mol_to_graph(mol2)

    # Generate subgraphs
    analyser = MolecularGraphAnalyser(calc.db)
    subgraphs1 = generate_all_subgraphs(g1)
    subgraphs2 = generate_all_subgraphs(g2)

    # Convert to hash sets for comparison
    hashes1 = Set(graph_to_canonical_hash(sg) for sg in subgraphs1)
    hashes2 = Set(graph_to_canonical_hash(sg) for sg in subgraphs2)

    # Calculate overlap metrics
    intersection = intersect(hashes1, hashes2)
    union_set = union(hashes1, hashes2)

    # Tanimoto coefficient (Jaccard similarity)
    tanimoto = length(intersection) / length(union_set)

    # Approximate GED based on subgraph overlap
    max_nodes = max(nv(g1), nv(g2))
    min_shared = length(intersection)
    max_possible_shared = min(length(hashes1), length(hashes2))

    # Normalised edit distance approximation
    if max_possible_shared == 0
        ged_approx = max_nodes
    else
        similarity_ratio = min_shared / max_possible_shared
        ged_approx = max_nodes * (1 - similarity_ratio)
    end

    # Calculate exact GED (simplified - would need more sophisticated algorithm)
    exact_ged = calculate_exact_ged(g1, g2)

    return (
        approximate_ged=round(ged_approx, digits=3),
        tanimoto_coefficient=round(tanimoto, digits=3),
        shared_subgraphs=length(intersection),
        total_subgraphs_mol1=length(hashes1),
        total_subgraphs_mol2=length(hashes2),
        exact_ged=exact_ged
    )
end

function calculate_exact_ged(g1::SimpleGraph, g2::SimpleGraph)
    """Calculate exact graph edit distance (simplified implementation)"""
    # This is a placeholder - exact GED is NP-hard
    # You might want to use specialized algorithms or approximations
    try
        # Simple heuristic based on node and edge differences
        node_diff = abs(nv(g1) - nv(g2))
        edge_diff = abs(ne(g1) - ne(g2))
        return node_diff + edge_diff
    catch e
        @warn "Could not calculate exact GED: $e"
        return nothing
    end
end

# Visualization functions (simplified)
function visualize_comparison(mol1::MoleculeInput, mol2::MoleculeInput, results)
    """Create a simple text-based comparison visualization"""
    println("\n" * "="^60)
    println("MOLECULAR COMPARISON RESULTS")
    println("="^60)
    println("Molecule 1: $(something(mol1.name, "Unknown")) ($(mol1.smiles))")
    println("Molecule 2: $(something(mol2.name, "Unknown")) ($(mol2.smiles))")
    println()
    println("Approximate GED: $(results.approximate_ged)")
    println("Exact GED: $(something(results.exact_ged, "N/A"))")
    println("Tanimoto Coefficient: $(results.tanimoto_coefficient)")
    println("Shared Subgraphs: $(results.shared_subgraphs)")
    println("Total subgraphs (mol1): $(results.total_subgraphs_mol1)")
    println("Total subgraphs (mol2): $(results.total_subgraphs_mol2)")
    println("="^60)
end

function plot_similarity_matrix(molecules::Vector{MoleculeInput},
    calculator::GraphEditDistanceCalculator)
    """Create a similarity matrix plot"""
    n = length(molecules)
    similarity_matrix = zeros(n, n)

    for i in 1:n
        for j in 1:n
            if i == j
                similarity_matrix[i, j] = 1.0
            elseif i < j
                results = calculate_ged_approximation(calculator, molecules[i], molecules[j])
                sim = results.tanimoto_coefficient
                similarity_matrix[i, j] = sim
                similarity_matrix[j, i] = sim
            end
        end
    end

    # Create heatmap
    mol_names = [something(mol.name, "Mol $i") for (i, mol) in enumerate(molecules)]

    p = heatmap(similarity_matrix,
        xticks=(1:n, mol_names),
        yticks=(1:n, mol_names),
        title="Molecular Similarity Matrix (Tanimoto Coefficients)",
        color=:viridis,
        aspect_ratio=:equal)

    # Add text annotations
    for i in 1:n, j in 1:n
        annotate!(j, i, text(string(round(similarity_matrix[i, j], digits=3)),
            :white, :center, 8))
    end

    return p
end

# Main execution
function main()
    """Main execution function with example usage"""

    # Initialize database and components
    db = SubgraphDatabase("molecular_subgraphs.db")
    analyser = MolecularGraphAnalyser(db)
    calculator = GraphEditDistanceCalculator(db)

    try
        # Example molecules (simplified SMILES for testing)
        examples = [
            MoleculeInput("CCCC", "Butane"),
            MoleculeInput("CCC", "Propane"),
            MoleculeInput("CC", "Ethane"),
            MoleculeInput("CCCCC", "Pentane"),
            MoleculeInput("CC(C)C", "Isobutane")
        ]

        # Populate database with subgraphs
        @info "Populating subgraph database..."
        for mol in examples
            try
                populate_database!(analyser, mol)
            catch e
                @warn "Failed to process molecule $(mol.name): $e"
            end
        end

        # Calculate graph edit distances between pairs
        @info "Calculating graph edit distances..."

        # Example comparison
        mol1 = examples[1]  # Butane
        mol2 = examples[5]  # Isobutane

        results = calculate_ged_approximation(calculator, mol1, mol2)
        visualize_comparison(mol1, mol2, results)

        # Additional comparisons
        println("\nAdditional Comparisons:")
        println("="^50)

        for i in 1:length(examples)
            for j in (i+1):length(examples)
                mol_a, mol_b = examples[i], examples[j]
                try
                    results = calculate_ged_approximation(calculator, mol_a, mol_b)
                    println("$(mol_a.name) vs $(mol_b.name): GED ≈ $(results.approximate_ged), " *
                            "Tanimoto = $(results.tanimoto_coefficient)")
                catch e
                    @warn "Failed comparison $(mol_a.name) vs $(mol_b.name): $e"
                end
            end
        end

        # Create similarity matrix plot
        try
            p = plot_similarity_matrix(examples, calculator)
            display(p)
        catch e
            @warn "Could not create similarity plot: $e"
        end

    catch e
        @error "Error in main execution: $e"
        rethrow(e)
    finally
        close(db)
    end
end

# Utility functions for interactive use
function quick_comparison(smiles1::String, smiles2::String,
    name1::String="Mol1", name2::String="Mol2")
    """Quick comparison of two molecules by SMILES"""
    db = SubgraphDatabase(":memory:")  # In-memory database
    calculator = GraphEditDistanceCalculator(db)

    try
        mol1 = MoleculeInput(smiles1, name1)
        mol2 = MoleculeInput(smiles2, name2)

        results = calculate_ged_approximation(calculator, mol1, mol2)
        visualize_comparison(mol1, mol2, results)

        return results
    finally
        close(db)
    end
end

# Export main functions
export MoleculeInput, SubgraphDatabase, MolecularGraphAnalyser
export GraphEditDistanceCalculator, main, quick_comparison

# Run main if this is the executed file
if abspath(PROGRAM_FILE) === @__FILE__
    main()
end