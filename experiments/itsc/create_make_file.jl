using Gadfly
using DataFrames

function p(n)
  if n < 10
    return "00$n"
  elseif n < 100
    return "0$n"
  end
  return "$n"
end

betas = collect(0:0.01:1.5)
index = 1
fd = open("Makefile", "w")
for beta in betas
  a = "all: exp$index\n"
  b = "exp$index:\n"
  c = "\t./experiments/itsc/bin/random_network -i 100000 -b $beta -o exp$(p(index))_$(beta).csv\n\n"
  write(fd, "$a$b$c")
  index = index + 1
end
close(fd)


betas = collect(0:0.01:1.5)
index = 1
for beta in betas
  fd = open("exp$(p(index))_$(beta).sh","w")
  a = "#!/bin/bash\n"
  b = "cd /usr/people/zahedi/experiments/random_network/\n"
  c = "/usr/people/zahedi/src/entropy/build/experiments/itsc/bin/random_network -i 1000 -b $beta -o exp$(p(index))_$(beta).csv\n\n"
  write(fd, "$a$b$c")
  index = index + 1
  close(fd)
end


files = readdir("/Users/zahedi/projects/builds/entropy-ninja/")
files = filter(x->contains(x,".csv"), files)

mean_values = zeros(length(files))
betas       = zeros(length(files))

for i = 1:length(files)
  a = readcsv("/Users/zahedi/projects/builds/entropy-ninja/$(files[i])")
  mean_values[i] = mean(a[2:end])
  betas[i] = a[1]
end

df = DataFrame()
df[:Beta] = betas
df[:Means] = mean_values

sort!(df, cols = [:Beta])

pp = plot(x=df[:Beta], y=df[:Means])
