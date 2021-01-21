struct PPR
    id::String
    motif_types::Vector{String}
    motifs::Vector{String}
    motif_scores::Vector{Float64}
end

struct TargetSequence
    id::String
    seq::String
    TargetSequence(i,s) = new(i,uppercase(s))
end

const PPR_motif_types = ["P","P1","L1","S1","SS","P2","L2","S2","E1","E2","Pi","P1i","L1i","S1i","SSi","P2i","L2i","S2i","E1i","E2i"]

const FORWARD_STRAND = 0x02
const REVERSE_STRAND = 0x01

function readPPRfinderMotifs(file)
    pprs = Vector{PPR}()
    current_id = ""
    current_motif_types = Vector{String}()
    current_motifs = Vector{String}()
    current_motif_scores = Vector{Float64}()
    for line in eachline(file)
       fields = split(line,"\t")
       if fields[1] ≠ current_id
           if !isempty(current_motifs)
                push!(pprs,PPR(current_id,current_motif_types,current_motifs,current_motif_scores))
                current_motif_types = Vector{String}()
                current_motifs = Vector{String}()
                current_motif_scores = Vector{Float64}()
           end
           current_id = fields[1]
       end
       if fields[9] in PPR_motif_types
           push!(current_motif_types,fields[9])
           push!(current_motifs,fields[7]*fields[8])
           push!(current_motif_scores,parse(Float64,fields[4]))
       end
    end
    if !isempty(current_motifs)
        push!(pprs,PPR(current_id,current_motif_types,current_motifs,current_motif_scores))
    end
    return pprs::Vector{PPR}
end

function extractFromMultiFasta(file,ids)
    seqs = Vector{TargetSequence}()
    current_id = nothing
    current_seq = IOBuffer()
    while !eof(file)
        line = readline(file)
        length(line) == 0 && continue
        if line[1] == '>'
            if !isnothing(current_id) && (isnothing(ids) || current_id in ids)
                push!(seqs,TargetSequence(current_id,String(take!(current_seq))))
            end
            current_id = split(line," ")[1][2:end]
        else
            print(current_seq, line)
        end
    end
    if !isnothing(current_id) && (isnothing(ids) || current_id in ids)
        push!(seqs,TargetSequence(current_id,String(take!(current_seq))))
    end
    return seqs
end

function readScoringTables(paths::Vector{String})
    scoring_tables = Dict{String,Dict{String,NamedTuple{(:A, :C, :G, :T, :U),Tuple{Float64,Float64,Float64,Float64,Float64}}}}()
    for path in paths
        if isdir(path) #check if file is directory, if so do this for all files in the directory
            files = readdir(path, join = true)
        else
            files = Vector{String}()
            push!(files,path)
        end
        for file in files
            scoring_table = Dict{String,NamedTuple{(:A, :C, :G, :T, :U),Tuple{Float64,Float64,Float64,Float64,Float64}}}()
            lines = readlines(file)
            motif_types = split(lines[1][18:end], r" |,|\t")
            for line in lines[3:end]
                fields = split(line,'\t')
                scoring_table[fields[1]] = (A=parse(Float64,fields[2]),C=parse(Float64,fields[3]),G=parse(Float64,fields[4]),T=parse(Float64,fields[5]),U=parse(Float64,fields[5]))
            end
            for motif_type in motif_types
                scoring_tables[motif_type] = scoring_table
            end
        end
    end
    return scoring_tables
end

const offsets = Dict("E2"=>2,"E1"=>3,"S2"=>4,"L2"=>5,"P2"=>6,"E2i"=>2,"E1i"=>3,"S2i"=>4,"L2i"=>5,"P2i"=>6)

function scoreAlignment(ppr,target,offset,scoring_tables)
    sum_score = 0.0
    for (i,(motif_type, motif)) in enumerate(zip(ppr.motif_types, ppr.motifs))
        st = get(scoring_tables, motif_type, nothing)
        isnothing(st) && continue
        combo = get(st,motif,nothing)
        isnothing(combo) && continue
        nt = target.seq[offset+i-1]
        score = getfield(combo,Symbol(nt))
        sum_score += score
    end
    return sum_score
end

function scoreAllAlignments(ppr,target,scoring_tables)
    maxoffset = length(target.seq)-length(ppr.motifs)+1
    scores = Vector{Float64}(undef,maxoffset)
    for offset in 1:maxoffset
        score = scoreAlignment(ppr,target,offset,scoring_tables)
        scores[offset] = score
    end
    return scores
end

const comp = Dict('G' => 'C', 'T' => 'A', 'A' => 'T', 'C' => 'G')

function revcomp(dna)
    return TargetSequence(dna.id*"_reverse",reverse(map(x -> comp[x], dna.seq)))
end

using Statistics
function scores2zscores(scores)
    mean = Statistics.mean(scores)
    std = Statistics.std(scores)
    zscores = zeros(Float64,length(scores))
    for (i,score) in enumerate(scores)
        zscores[i] = (score-mean)/std
    end
    return mean,std,zscores
end

using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--reference", "-r"
            help = "path to reference sequence for calculating score distribution"
            default = nothing
        "--edit_site", "-e"
            help = "edit site position within the sequence(s)"
            arg_type = Int
            default = 0
        "--strand", "-s"
            help = "which target strand(s) to match to"
            arg_type = String
            default = "F"
        "--outfile", "-o"
            help = "path to outfile for saving results"
            arg_type = String
            default = nothing
        "PPR(s)"
            help = "path to file of PPRfinder motifs to search with"
            arg_type = String
            required = true
        "target(s)"
            help = "path to file of FASTA format target sequence(s) to search in"
            arg_type = String
            required = true
        "scoring_table(s)"
            help = "path to scoring table or directory of scoring table files"
            nargs = '+'
            arg_type = String
            required = true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()

    local pprs
    open(parsed_args["PPR(s)"]) do infile_motifs
        pprs = readPPRfinderMotifs(infile_motifs)
    end

    scoring_tables = readScoringTables(parsed_args["scoring_table(s)"])

    using_reference = false
    ref_file = parsed_args["reference"]
    if !isnothing(ref_file)
        using_reference = true
        ref_means = Vector{Float64}(undef, length(pprs))
        ref_stds = Vector{Float64}(undef, length(pprs))
        open(ref_file) do reff
            reference = extractFromMultiFasta(reff,nothing)[1]
            for (i, ppr) in enumerate(pprs)
                fscores = scoreAllAlignments(ppr,reference,scoring_tables)
                rscores = scoreAllAlignments(ppr,revcomp(reference),scoring_tables)
                mean,std,zscores = scores2zscores(vcat(fscores,rscores))
                ref_means[i] = mean
                ref_stds[i] = std
            end
        end
    end

    local targets
    open(parsed_args["target(s)"]) do infile_targets
        targets = extractFromMultiFasta(infile_targets,nothing)
    end

    strandarg = parsed_args["strand"]
    strand = 0x00
    if length(strandarg) > 3
        if occursin(r"forward"i, strandarg); strand |= 0x02; end
        if occursin(r"reverse"i, strandarg); strand |= 0x01; end
        if occursin(r"both"i, strandarg); strand |= 0x03; end
    else
        if occursin(r"f|\+"i, strandarg); strand |= 0x02; end
        if occursin(r"r|-"i, strandarg); strand |= 0x01; end
    end

    if strand & REVERSE_STRAND ≠ 0
        rev_targets = Vector{TargetSequence}()
        for target in targets
            push!(rev_targets, revcomp(target))
        end
        if strand & FORWARD_STRAND ≠ 0
            targets = append!(targets,rev_targets)
        else
            targets = rev_targets
        end
    end

    #create outfile
    outfile = parsed_args["outfile"]
    if !isnothing(outfile)
        io = open(outfile, "w")
    else
        io = stdout
    end

    edit_site = parsed_args["edit_site"]

    write(io, join(["protein_ID","target_ID","start","end","aligned_nucleotides","score"],"\t"))
    if using_reference; write(io, "\tZscore"); end
    write(io, "\n")

    for (i, ppr) in enumerate(pprs), target in targets
        if edit_site > 0
            offset = edit_site - offsets[last(ppr.motif_types)] - length(ppr.motifs) + 1
            score = scoreAlignment(ppr, target, offset, scoring_tables)
            print(io, join([ppr.id, target.id, offset, offset+length(ppr.motifs)-1, target.seq[offset:offset+length(ppr.motifs)-1], score], "\t"))
            if using_reference
                print(io, "\t", (score-ref_means[i])/ref_stds[i])
            end
            println(io)
        else
            scores = scoreAllAlignments(ppr,target,scoring_tables)
            for (pos,score) in enumerate(scores)
                print(io, join([ppr.id, target.id, pos, pos+length(ppr.motifs)-1, target.seq[pos:pos+length(ppr.motifs)-1], score], "\t"))
                if using_reference
                    print(io, "\t", (score-ref_means[i])/ref_stds[i])
                end
                println(io)
            end
        end
    end
    close(io)
end

main()
