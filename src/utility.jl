function warn_ifdefined_quote(s::Symbol)

    warnstr = """
    Symbol $s is already defined in current scope! It might be overwritten now.
    """
    return quote
        if @isdefined $s
            @warn $warnstr
        end
    end    
end

