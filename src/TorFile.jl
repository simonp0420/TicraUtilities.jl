const LP = '('
const RP = ')'


"""
    TorObj

Struct containing a single, general, TOR object.

### Fields
* `name::String`: The name of the TOR object. Example: "my_rim".
* `objtype::String`: The type of the TOR object. Example: "tabulated_rim_xy".
* `propname::Vector{String}`: Names of object properties. Example: 
  ["file_name", "unit", "number_of_points", "translation", "polar_origin"]
* `propval::Vector{String}`: Values of object properties corresponding to the names in `propnames`. Example:
  ["h_9m_scalloped_rim.xyz", "in", "112", "struct(x: 0.0 in, y: 0.0 in)", 
  "struct(status: automatic, x: 0.0 in, y: 0.0 in)"]

Note that the "values" in `propval` are simply strings that must be parsed for any numeric values
that might be present.
"""
struct TorObj
    name::String  # Name of the object
    objtype::String  # Type of the object
    propname::Vector{String}  # names of object properties
    propval::Vector{String}  # values of object properties
end

"""
    get_name(obj::TorObj) -> name::String

Return the object name.
"""
get_name(obj::TorObj) = obj.name

"""
    get_objtype(obj::TorObj) -> objtype::String

Return the Ticra object type.
"""
get_objtype(obj::TorObj) = obj.objtype

"""
    get_propname(obj::TorObj) -> propname::Vector{String}

Return the vector of object property names.
"""
get_propname(obj::TorObj) = obj.propname

"""
    get_propval(obj::TorObj) -> propval::Vector{String}

Return the vector of object property values.
"""
get_propval(obj::TorObj) = obj.propval


Base.show(io::IO, mime::MIME"text/plain", obj::TorObj) = _write_tor_object(io, obj)
Base.show(io::IO, obj::TorObj) = print(io, obj.name, "::", obj.objtype)

"""
    write_torfile(fname::AbstractString, obj::TorObj; modestr = "w")
    write_torfile(fname::AbstractString, objs::Vector{TorObj}; modestr = "w")

Write one or multiple `TorObj` objects to a Ticra-compatible TOR (Ticra Object Repository) file.
With the default value of `modestr = "w"`, the file is created if it doesn't already exist, or, if it 
already exists, it is truncated to zero length before writing the TOR object(s).
By setting the value of `modestr` to `"a"` one can append the object(s) to an existing file.
"""
function write_torfile(fname::AbstractString, objs::Vector{TorObj}; modestr = "w")
    open(fname, modestr) do io
        for obj in objs
            _write_tor_object(io, obj)
        end
    end
    return nothing
end

write_torfile(fname::AbstractString, obj::TorObj; kwargs...) = write_torfile(fname, [obj]; kwargs...)

"""
    _write_tor_object(io::IO, obj::TorObj)

Write a tor object to a previously opened tor output file.
"""
function _write_tor_object(iu::IO, obj::TorObj)

    println(iu, obj.name, "   ", obj.objtype, " (")

    for (ip, name) in enumerate(obj.propname)
        print(iu, "   ", name, " : ")
        L = length(name) + 6  # Length of string written so far
        line = obj.propval[ip] # Remainder to be written
        while length(line) > 0
            if L + length(line) ≤ 80
                print(iu, line) # Write the whole thing
                break
            else
                k = findlast(view(line, 1:80-L), ",")
                if isnothing(k)
                    print(iu, line)
                    break
                end
                print(iu, view(line, 1:k))
                line = strip(line[k+1:length(line)])
                if length(line) > 0
                    println(iu, "")
                    print(iu, " "^L)
                end
            end
        end
        ip < length(obj.propname) && println(iu, ",")
    end
    println(iu, ")")
    return nothing
end

"""
    contents::Vector{String} = _read_tor_file(filename::AbstractString)

Read an entire TOR file into memory.  Each line of the file is stored in
an element of contents.  Comments are removed, as is white space at start 
and end of lines.  Tabs are replaced by blanks.
"""
function _read_tor_file(filename)
    CE = "*/" # comment end
    CS = "/*" # comment start
    in_block_comment = false  # /* */-style comment flag

    content::Vector{String} = open(readlines, filename, "r")


    for (k, line) in enumerate(content)
        line = strip(replace(line, '\t' => ' '))
        # Strip off material after // comment marker:
        ir = findfirst("//", line)
        isnothing(ir) || (line = strip(line[1:(first(ir)-1)]))

        if in_block_comment
            ir = findfirst(CE, line) # Test if comment ends on this line:
            if isnothing(ir)
                line = ""  # Erase this entire line since entire line is within comment region
            else
                line = strip(line[(last(ir)+1):end]) # Throw away commented portion
                in_comment = false  # toggle comment state flag
            end
        else
            ir = findfirst(CS, line) # Test if comment starts on this line
            if !isnothing(ir)
                # This line has a comment start '/*' in it.  Check if the 
                # comment terminates on the same line:
                i1r = findfirst(CE, line)
                if !isnothing(i1r)
                    @assert i1r[1] > ir[end] "*/ precedes /* on line $(k) of tor file "
                    # Comment terminates on this line.
                    line = line[1:first(ir)-1] * line[last(i1r)+1:end]  # Remove the comment
                else
                    # Comment does not terminate on this line.
                    line = line[1:first(ir)-1]  # Remove the comment portion
                    in_comment = true  # toggle comment state
                end
            end
        end
        content[k] = strip(line)
    end
    return content
end

"""
    (token, i1) = _next_token_from_i1(line, i1)

Alters i1 to point to next character after end of token.

"""
function _next_token_from_i1(line, i1)
    while i1 ≤ length(line) && isspace(line[i1])
        i1 += 1
    end
    i2 = findnext(x -> x ∈ (' ', ':'), line, i1) - 1
    if isnothing(i2)
        return (line[i1:end], length(line))
    else
        return (line[i1:i2], i2 + 1)
    end
end

"""
    (i1, propnames, propvals) = _parse_properties!(line, i1)

On entry, `i1` is pointing to the opening parenthesis of the property block.  
On exit, `i1` is pointing to next character after end (i.e., closing parenthesis)
of properties block.
"""
function _parse_properties(line, i1)
    propnames = String[]
    propvals = String[]
    i1 += 1
    while true
        while isspace(line[i1])
            i1 += 1
        end
        if line[i1] == RP
            i1 += 1
            return (i1, propnames, propvals)
        end
        (propname, i1) = _next_token_from_i1(line, i1)
        push!(propnames, propname)
        while isspace(line[i1])
            i1 += 1
        end
        @assert line[i1] == ':'
        i1 += 1
        i2 = _find_next_comma_or_level0_rp(line, i1)
        push!(propvals, strip(line[i1:i2-1]))
        i1 = i2 + 1
        if line[i2] == RP
            return (i1, propnames, propvals)
        end
    end
end

"""
    i2 = _find_next_comma_or_level0_rp(line, i1)

Find next comma or right parenthesis not preceded by a left parenthesis
"""
function _find_next_comma_or_level0_rp(line, i1)
    plevel = 1 # paren level
    i2 = i1 + 1
    while i2 < length(line)
        i2 += 1
        if line[i2] == LP
            plevel += 1
        elseif line[i2] == RP
            plevel -= 1
            plevel < 1 && (return i2)
        end
        if plevel == 1 && line[i2] == ','
            return i2
        end
    end
end




"""
    read_torfile(torfile::AbstractString)

Return a vector of TorObj objects found in a TOR file
"""
function read_torfile(torfile::AbstractString)
    i1 = 1
    line = strip(join(_read_tor_file(torfile), ' ')) # Create a single, long string to process
    objects = TorObj[]
    while i1 < length(line)
        # i1 is pointing to the first character of the object name
        (name, i1) = _next_token_from_i1(line, i1)
        isempty(name) && break
        (tortype, i1) = _next_token_from_i1(line, i1)
        @assert !isempty(tortype) "Empty tor type following object name!"
        push!(objects, TorObj(name, tortype, String[], String[]))
        i1 = findnext(LP, line, i1)  # Start of property list
        isnothing(i1) && (return objects)
        (i1, propnames, propvals) = _parse_properties(line, i1)
        append!(objects[end].propname, propnames)
        append!(objects[end].propval, propvals)
    end
    return objects
end

"""
    i = search_name_index(objs::Vector{TorObj}, name::AbstractString)

Return index i of object in objs with objs[i].name == name
"""
function search_name_index(objs::Vector{TorObj}, name::AbstractString)
    findfirst(==(name), (o.name for o in objs))
end

"""
    (x, y, xunit, yunit) = parse_tor_xy_struct(s::AbstractString)

The input string is an entry in a TorObj propval field, consisting
of, e.g., "struct(x: 75.9880 in, y: 0.0 in)".  This function
parses out and returns the floating point x and y values and their 
respective unit strings.
"""
function parse_tor_xy_struct(s::AbstractString)
    cf = r"([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"  # This RE captures a float
    (x, y) = (parse(Float64, s[rng]) for rng in findall(cf, s))
    s2 = split(s, ',')
    indx = findfirst(str -> occursin("x:", str), s2)
    isempty(indx) && error("Unable to find \"x:\" in string")
    xunit = split(s2[indx])[end]
    indx = findfirst(str -> occursin("y:", str), s2)
    isempty(indx) && error("Unable to find \"y:\" in string")
    yunit = split(s2[indx])[end]
    yunit[end] == ')' && (yunit = yunit[1:end-1])
    return (x, y, xunit, yunit)
end

