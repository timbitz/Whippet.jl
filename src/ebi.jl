

type EBIResponse
   fastq_1_url::String
   fastq_2_url::String
   paired::Bool
   success::Bool
end

function ident_to_fastq_url( ebi_id::String )
   request_url = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=$ebi_id&result=read_run&fields=fastq_ftp&display=txt"
   res = HTTP.get(request_url)
   res.status == 200 || error("Couldn't retrieve information for $ebi_id from ebi.ac.uk!")
   buf = res.body
   header = chomp(readline(buf))
   if eof(buf)
      error("No fastq data is associated with $ebi_id according to ebi.ac.uk!")
   end
   urls   = split(chomp(readline(buf)), ';', keep=false)
   if length(urls) >= 2
      success = res.status == 200 && urls[1] != "" && urls[2] != "" ? true : false
      ret = EBIResponse(string(urls[1]), string(urls[2]), true, success)
   else
      success = res.status == 200 && urls[1] != "" ? true : false
      ret = EBIResponse(string(urls[1]), "", false, success)
   end
   ret
end
