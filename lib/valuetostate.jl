module ValueToState

function getState(value, downregulatedCutoff, upregulatedCutoff)
    if upregulatedCutoff < value
        return "Up Regulated"
    elseif downregulatedCutoff > value
        return "Down Regulated"
    else
        return "Normal"
    end
end

end
