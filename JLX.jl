import Base: BitUnsigned

mutable struct JLX{T<:BitUnsigned}
    fork_lcg::T
    main_lcg::T
    main_xsh::Tuple{T,T}
end

JLX{T}() where {T<:BitUnsigned} =
    JLX{T}(rand(T), rand(T), (rand(T), rand(T)))
JLX(T::Type{<:BitUnsigned} = UInt64) = JLX{T}()

function next(rng::JLX)
    lcg = rng.main_lcg
    xsh = rng.main_xsh
    out = @inline mix(xsh[1], lcg)
    xsh = @inline update_xsh(xsh)
    lcg = @inline update_lcg(lcg)
    rng.main_xsh = xsh
    rng.main_lcg = lcg
    return out
end

# random additive offsets
const A1 = 0x7d798fbd7805bf42
const A2 = 0xda7e49f8f5af63f6
const A3 = 0xe054ee0b61f0891d
const A4 = 0x39455e68c9093594

function fork(rng::JLX{T}) where {T<:BitUnsigned}
    flcg = rng.fork_lcg
    flcg = @inline update_lcg(flcg)
    rng.fork_lcg = flcg
    xsh1 = @inline mix(flcg + (A1 % T), rng.main_xsh[1] + (A2 % T))
    xsh2 = @inline mix(flcg + (A3 % T), rng.main_xsh[2] + (A4 % T))
    return JLX(flcg, rng.main_lcg, (xsh1, xsh2))
end

# xoshiro/xorshift implementations

update_xsh(xsh::Tuple) = update_xsh(xsh...)

# xoroshiro128
function update_xsh(s0, s1)
	s1 ⊻= s0
	s0 = bitrotate(s0, 24) ⊻ s1 ⊻ (s1 << 16)
	s1 = bitrotate(s1, 37)
    s0, s1
end

# LCG implementations

update_lcg(lcg::UInt64) = lcg*0xd1342543de82ef95 + 1

# non-linear word mixing

function rxs_shift(a::BitUnsigned)
    b = sizeof(a) << 3
    c = b - leading_zeros(b) - 2
    a >> (b - c) + c
end

# 2/3 of the bits
xs_shift(a::BitUnsigned) = div(sizeof(a) << 4, 3, RoundNearest)

function mix(a::T, b::T) where {T<:BitUnsigned}
    c = a ⊻ (b >> rxs_shift(a))
    d = b ⊻ (a >> rxs_shift(b))
    e = ((c*d) << 1) + c + d
    e ⊻= e >> xs_shift(e)
end

function mix(a::UInt64, b::UInt64)
    c = ((a*b) << 1) + a + b
    c ⊻= c >> 43
end
