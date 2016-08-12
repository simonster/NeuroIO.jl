module Struct
export @struct
import Compat.String

function struct_string(bytes::Vector{UInt8})
    last_byte = findfirst(bytes, 0)
    String(last_byte == 0 ? bytes : bytes[1:last_byte-1])
end

macro struct(typename, contents)
    if contents.head != :block
        error("Invalid struct declaration")
    end

    read_blk = Expr(:call, esc(typename))
    write_blk = Expr(:block)

    for typedecl in contents.args
        if typedecl.head == :line
            continue
        elseif typedecl.head != :(::)
            error("invalid struct declaration")
        end
        fieldname = typedecl.args[1]
        fieldtype = typedecl.args[2]
        field = :(x.$(fieldname))
        if isa(fieldtype, Expr) && fieldtype.head != :curly
            if fieldtype.head != :call
                error("invalid struct declaration")
            end
            # Type has size parameters

            if fieldtype.args[1] in (:ASCIIString, :UTF8String, :String)
                typedecl.args[2] = fieldtype.args[1]
                push!(read_blk.args, :(struct_string(read(io, UInt8, ($(fieldtype.args[2:end]...))))))
                push!(write_blk.args, quote
                    if sizeof($field) > $(fieldtype.args[2])
                        throw(ArgumentError($("size of $(fieldtype.args[1]) must be less than $(fieldtype.args[2])")))
                    end
                    write(io, $field)
                    for i = sizeof($field)+1:$(fieldtype.args[2])
                        write(io, UInt8(0))
                    end
                end)
            else
                typedecl.args[2] = :(Array{$(fieldtype.args[1]), $(length(fieldtype.args)-1)})
                push!(read_blk.args, :(read(io, $(esc(fieldtype.args[1])), ($(fieldtype.args[2:end]...)))))
                push!(write_blk.args, quote
                    if size($field) != $(tuple(fieldtype.args[2:end]...))
                        throw(ArgumentError($("size of $(fieldname) must be $((fieldtype.args[2:end]...,)), but was ")*string(size($field))))
                    end
                    write(io, $field)
                end)
            end
        else
            push!(read_blk.args, :(read(io, $fieldtype)))
            push!(write_blk.args, :(write(io, $field)))
        end
    end
    quote
        type $(esc(typename))
            $(esc(contents))
        end
        Base.read(io::IO, ::Type{$(esc(typename))}) = $read_blk
        Base.write(io::IO, x::$(esc(typename))) = $write_blk
    end
end
end
