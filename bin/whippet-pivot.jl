#!/usr/bin/env julia
# Tim Sterne-Weiler 2021

using Pkg

const dir = abspath( splitdir(@__FILE__)[1] )
const ver = chomp(readline(open(dir * "/VERSION")))

start = time_ns()
println( stderr, "Whippet $ver loading... " )

Pkg.activate(dir * "/..")

using Whippet
using Random
using ArgParse
using Glob
using BufferedStreams
using Libz
using CSV
using DataFrames


function parse_cmd()
  s = ArgParseSettings()
  @add_arg_table s begin
    "--files", "-f"
      help     = "Pattern to glob.psi[.gz] (common-filename-segment [*.psi*]), or comma delimited list of filenames. ie. (-a sample_a) would work for sample_a-rep1.psi.gz,sample_a-rep2.psi.gz,..."
      arg_type = String
      default  = ""
#      required = true
    "--directory", "-d"
      help     = "Directory to search for file patterns or list in -a and -b"
      arg_type = String
      default  = "."
    "--gene-centric"
    	help		= "Pivot sample_centric .psi[.gz] files into gene_centric .gpsi.gz files"
    	action   = :store_true
    "--output-extension"
      help     = "Default output extension"
      arg_type = String
      default  = ".gpsi.gz"
    "--matrix"
      help     = "Make a matrix of nodes vs. samples for a particular value type (psi, entropy)"
      arg_type = String
    "--bufsize"
      help     = "Buffer size"
      arg_type = Int64
      default  = 1_048_576 #1_024*8
  end
  return parse_args(s)
end

function retrievefilelist( pattern::String, dir::String )
   list = Vector{String}()
   if something(findfirst(isequal(','), pattern), 0) > 0
      tmp = split( pattern, ',', keepempty=false )
   else
      pat = length(pattern) > 0 ? "*" * pattern : pattern
      tmp = glob( pat * "*psi*", dir )
   end
   # now clean the return
   for file in tmp
      push!( list, string(file) )
   end
   list
end

readlinesplit( stream ) =  split( readline(stream), '\t' )

function gene_centric( streams::Vector{BufferedStreams.BufferedInputStream},
                       names::Vector{String};
                       var_column::Int64=0,
                       output_ext::String=".gpsi.gz" )
   
   curline  = Vector{Vector{SubString{String}}}(undef, length(streams))
   metadata = Vector{String}()
   header   = ""

   # initialize header lines
   for (i,s) in enumerate(streams)
   	header = readline(s)
   	curline[i] = readlinesplit(s) # first psi line
   end
   curgene = curline[1][1]
   ostream = ZlibDeflateOutputStream(open(curgene * output_ext, "w"))
   write(ostream, "Sample\t" * header * "\n")
   println(stderr, "Pivoting into " * curgene * output_ext)

   while !eof(streams[1])
   	for (i,s) in enumerate(streams)
      	while curline[i][1] == curgene && !eof(s)
      		tab_write(ostream, names[i])
      		tab_write(ostream, curline[i])

      		curline[i] = readlinesplit(s)
      	end
      end
      close(ostream)
      curgene = curline[1][1]
   	  ostream = ZlibDeflateOutputStream(open(curgene * output_ext, "w"))
      write(ostream, "Sample\t" * header * "\n")
      println(stderr, "Pivoting into " * curgene * output_ext)
   end
end


function psi_to_dataframe( files::Vector{String}, summarize_type=true )
   vdf = Vector{DataFrame}()

   for i in 1:length(files)
      if summarize_type
         full = CSV.read(bufferedinput(files[i]), DataFrame, delim='\t', types=Dict(:Sample=>CategoricalArray, :Gene=>CategoricalArray)) |>
                             x->select(x, [:Sample, :Gene, :Node, :Type, :Psi, :Entropy]) |>
                             x->filter(row -> row.Psi != "NA", x) |>
                             x->filter(row -> row.Entropy != "NA", x)
         full[!,:Psi ]     = tryparse.(Float64,full[:,:Psi ])
         full[!,:Entropy ] = tryparse.(Float64,full[:,:Entropy ])
         full[!,:Sample ] = CategoricalArray(full[!,:Sample ])
         full[!,:Gene ]   = CategoricalArray(full[!,:Gene ])
         full[!,:Type ]   = CategoricalArray(full[!,:Type ])
         reg_df = filter(row -> row.Type in ["CE", "AA", "AD", "RI"], full)
         #reg_df[!,:Node ]  = parse.(Int64,reg_df[:,:Node])
         #sum_df = filter(row -> row.Type in ["XRI", "XDI", "XAI", "XDE", "XAE", "XEI", "XCE"], full) |>
         #                    x->groupby(x, [:Sample, :Gene, :Type]) |>
         #                    x->combine(x, [:Psi, :Entropy] .=> mean)
         push!(vdf, reg_df)
         #push!(vdf, sum_df)
      else
         temp = CSV.read(bufferedinput(files[i]), DataFrame, delim='\t', types=Dict(:Psi => Float64, :Entropy => Float64)) |>
                    x->select(x, [:Sample, :Gene, :Node, :Type, :Psi, :Entropy]) |>
                    x->filter(row -> row.Psi != "NA", x) |>
                    x->filter(row -> row.Entropy != "NA", x)
         temp[!,:Psi ] = tryparse.(Float64,temp[:,:Psi ])
         temp[!,:Entropy ] = tryparse.(Float64,temp[:,:Entropy ])
         push!(vdf, temp)
      end
   end
   vcat( vdf... )
end

unstack_psi( df::DataFrame )     = unstack_col( df, :Psi )
unstack_entropy( df::DataFrame ) = unstack_col( df, :Entropy )

function unstack_col( df::DataFrame, feature::Symbol=:Psi )
   df  = transform(df, [:Gene, :Node] => ByRow((Gene, Node) -> string(Gene, "_", Node)) => :GeneNode) |>
         x->select(x, [:Sample, :GeneNode, feature])

   udf = unstack(df, :Sample, :GeneNode, feature, allowduplicates=true)
   udf
end

function summarize_type( df::DataFrame )
   gdf = groupby(df, [:Sample, :Gene, :Node, :Type]) |>
      x->combine(gdf, :Psi => mean)

end

function main()
   args  = parse_cmd()
   println(stderr, " $( round( (time_ns()-start)/1e9, digits=6 ) ) seconds" )

   dir   = fixpath( args["directory"] )
   files = retrievefilelist( args["files"], dir )
	fnames = map(x->last(splitdir(x)), files) |>
	      x->map(y->String(first(split(y, r".(psi|isoform.tpm|gene.tpm)"))), x)

   println(stderr, "Loading files to pivot: $(join(map(basename, files), '\n'))")

   if length(files) <= 0
      error("Unable to find any files! length(files) == $(length(files))")
   end

   fstreams = open_streams( files, bufsize=args["bufsize"] )

   if args["gene-centric"]
   	println(stderr, "Pivoting to gene-centric psi files...")
   	gene_centric(fstreams, fnames)
   end

   println(stderr, "Whippet $ver done." )
end

@timer main()
