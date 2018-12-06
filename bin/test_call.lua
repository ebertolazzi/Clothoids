-- print contents of a table, with keys sorted. second parameter is optional, used for indenting subtables
function dump(t,indent)
    local names = {}
    if not indent then indent = "" end
    for n,g in pairs(t) do
        table.insert(names,n)
    end
    table.sort(names)
    for i,n in pairs(names) do
        local v = t[n]
        if type(v) == "table" then
            if(v==t) then -- prevent endless loop if table contains reference to itself
                print(indent..tostring(n)..": <-")
            else
                print(indent..tostring(n)..":")
                dump(v,indent.."   ")
            end
        else
            if type(v) == "function" then
                print(indent..tostring(n).."()")
            else
                print(indent..tostring(n)..": "..tostring(v))
            end
        end
    end
end

function pippo(a)
  dump(a,"");
  res = {
    sol = a[1]*a[2],
    cpy = a[4],
    aaa = {
      vec = {1,2,3,4},
      map = { a = 1, b = 2, c = 3 }
    },
  }
  return res;
end