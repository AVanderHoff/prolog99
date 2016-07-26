
len([],0).
len([_|T],N):-len(T,X),N is X+1.
myreverse(X,Y):-rev(X,Y,[]).
rev([],L,L):-!.
rev([H|T],L,A):-rev(T,L,[H|A]).
car([X|_],X).
my_last(X,[X]).
my_last(X,[_|L]):-my_last(X,L).
last_but_one(X,[X,_]).
last_but_one(X,[_,Y|Ys]):-last_but_one(X,[Y|Ys]).

elementat(X,[X|_],1).
elementat(X,[_|L],K):-K >1,K1 is K - 1,elementat(X,L,K1).

length(X,[],X).
length(X,[_|T],K):- K1 is K + 1,length(X,T,K1).

theirlength([],0).
theirlength([_|L],N):-theirlength(L,N1),N is N1 + 1.
eq([],[]).
eq([X|Tx],[Y|Ty]):- X = Y ,eq(Tx,Ty).

palindrome(X):- reverse(X,X).

my_flatten(X,[X]):- \+ is_list(X).
my_flatten([],[]).
my_flatten([X|Xs],Zs):-my_flatten(X,Y),my_flatten(Xs,Ys),append(Y,Ys,Zs).

compress([],Y,X):-reverse(Y,X).
compress([H|T],Acc,Z):-member(H,T),compress(T,Acc,Z).
compress([H|T],Acc,Z):-not(member(H,T)),compress(T,[H|Acc],Z).

theircompress([],[]).
theircompress([X],[X]).
theircompress([X,X|Xs],Zs):-theircompress([X|Xs],Zs).
theircompress([X,Y|Ys],[X|Zs]):- X\=Y ,theircompress([Y|Ys],Zs).


pack([],[]).
pack([X|Xs],[Z|Zs]):-transfer(X,Xs,Ys,Z),pack(Ys,Zs).


transfer(X,[],[],[X]).
transfer(X,[Y|Ys],[Y|Ys],[X]):-X\=Y.
transfer(X,[X|Xs],Ys,[X|Zs]):-transfer(X,Xs,Ys,Zs).

encode(X,Y):-pack(X,Z),packlists(Z,[],Y).
packlists([],S,W):-reverse(S,W).
packlists([H|T],S,W):-car(H,A),length(H,Q),Q==1,packlists(T,[A|S],W).
packlists([H|T],S,W):-car(H,A),length(H,Q),packlists(T,[[Q,A]|S],W).

theirEncode(L1,L2):-pack(L1,L),transform(L,L2).
transform([],[]).
transform([[X|Xs]|Ys],[[N,X]|Zs]):-length([X|Xs],N),transform(Ys,Zs).

theirModEncode(L1,L2):-theirEncode(L1,L),strip(L,L2).
strip([],[]).
strip([[1,X]|Ys],[X|Zs]):-strip(Ys,Zs).
strip([[N,X]|Ys],[[N,X]|Zs]):-N > 1, strip(Ys,Zs).

decode([],[]).
decode([X|Ys],[X|Zs]):- \+ is_list(X),decode(Ys,Zs).
decode([[1,X]|Ys],[X|Zs]):-decode(Ys,Zs).
decode([[N,X]|Ys],[X|Zs]):-N > 1, N1 is N - 1, decode([[N1,X]|Ys],Zs).

encode_direct([],[]).
encode_direct([X|Xs],[Z|Zs]):-count(X,Xs,Ys,1,Z),encode_direct(Ys,Zs).

count(X,[],[],1,X).
count(X,[],[],N,[N,X]):-N > 1.
count(X,[Y|Ys],[Y|Ys],N,[N,X]):-N > 1, X \= Y.
count(X,[X|Xs],Ys,K,T):- K1 is K + 1, count(X,Xs,Ys,K1,T).

dupli([],[]).
dupli([X|Tx],[X,X|Ty]):-dupli(Tx,Ty).
dupli_R(L1,N,L2):-dupli_R(L1,N,L2,N).
dupli_R([],_,[],_).
dupli_R([_|Xs],N,Ys,0):-dupli_R(Xs,N,Ys,N).
dupli_R([X|Xs],N,[X|Ys],K):- K > 0, K1 is K - 1 ,dupli_R([X|Xs],N,Ys,K1).

drop(L1,N,L2):-drop(L1,N,L2,N).
drop([],_,[],_).
drop([_|Xs],N,Ys,0):-drop(Xs,N,Ys,N).
drop([X|Xs],N,[X|Ys],K):- K > 0, K1 is K - 1, drop(Xs,N,Ys,K1).



split(Xs,0,[],Xs).
split([X|Xs],N,[X|Ys],Z):- N > 0, N1 is N - 1, split(Xs,N1,Ys,Z).

slice(L1,N,N1,L2):-slicehelp(L1,N,N1,1,L2).
slicehelp(_,_,W,K,[]):- K is W + 1.
slicehelp([_|Xs],N,N1,W,Ys):- W < N , W1 is W + 1, slicehelp(Xs,N,N1,W1,Ys).
slicehelp([X|Xs],N,N1,W,[X|Ys]):- W1 is W + 1, slicehelp(Xs,N,N1,W1,Ys).

theirslice([X|_],1,1,[X]).
theirslice([X|Xs],1,K,[X|Ys]):- K > 1,K1 is K - 1, theirslice(Xs,1,K1,Ys).
theirslice([_|Xs],I,K,Ys):- I > 1, I1 is I - 1, K1 is K - 1,
theirslice(Xs,I1,K1,Ys).

rotate(L2,N,L3):- N > 0, split(L2,N,X,Y), append(Y,X,L3).
rotate(L2,N,L3):- N < 0, length(L2,W), Z is N * -1, A is W - Z,
	split(L2,A,X,Y), append(Y,X,L3).

theirrotate([],_,[]):-!.
theirrotate(L1,N,L2):- length(L1,NL1),N1 is N mod NL1, split(L1,N1,S1,S2),
	append(S2,S1,L2).

removefirst([_|Tx],Tx).
first([X|_],X).
remove(E,L2,N,L3):- N1 is N - 1, split(L2,N1,Z,W),first(W,E), removefirst(W,Y),   append(Z,Y,L3).

%theirremove(X,[_|Xs],1,Xs).
% theirremove(X,[Y|Xs],K,[Y|Ys]):- K > 1, K1 is K - 1,
% theirremove(X,Xs,K1,Ys).

insert(X,L1,N,L2):-N1 is N - 1, split(L1,N1,Z,W) ,append(Z,[X|W],L2).
theirinsert(X,L,K,R):-theirremove(X,R,K,L).
range(N,W,[]):-	L is W + 1, N=L.
range(N,N1,[N|L]):- N2 is N + 1, range(N2,N1,L).

theirrange(I,I,[I]).
theirrange(I,K,[I|L]):- I < K, I1 is I + 1, range(I1,K,L).

rnd_select(_,0,[]).
rnd_select(X,N,[W|L]):- N1 is N - 1, length(X,Z),Q is Z - 1,random(0,Q,Y),
	nth0(Y,X,W), rnd_select(X,N1,L).

rnd_select2(X,Y,L):- range(1,Y,W),rnd_select(W,X,L).


rnd_perm(X,L):-length(X,Z),rnd_select(X,Z,L).

combination(0,_,[]).
combination(K,L,[X|Xs]):- K > 0 , el(X,L,R), K1 is K - 1, combination(K1,R,Xs).
el(X,[X|L],L).
el(X,[_|L],R):-el(X,L,R).

group3(X,G1,G2,G3):- combination(2,X,G1),subtract(X,G1,Y),combination(3,Y,G2),
	subtract(Y,G2,G3).
split3([X,Y,Z|_],X,Y,Z).
group32(X,Z,[G1,G2|[G3]]):-split3(Z,A,B,_),combination(A,X,G1),subtract(X,G1,Y),
	combination(B,Y,G2),subtract(Y,G2,G3).

selectN(0,_,[]):-!.
selectN(N,L,[X|S]):-N > 0, el(X,L,R),N1 is N - 1, selectN(N1,R,S).

lsort(In,Out,Dir):-add_key(In,KL,Dir),keysort(KL,SK),rem_key(SK,Out).
add_key([],[],_).
add_key([X|Xs],[L-p(X)|Ys],asc):-!,length(X,L),add_key(Xs,Ys,asc).
add_key([X|Xs],[L-p(X)|Ys],desc):-length(X,L1), L is -L1, add_key(Xs,Ys,desc).
rem_key([],[]).
rem_key([_-p(X)|Xs],[X|Ys]):-rem_key(Xs,Ys).


prime(Y):- primechecker(2,Y,Y*Y).
primechecker(X,Y,Z):- X >= Z ; S is mod(X,Y), S == 0.
primechecker(X,Y,Z):- X < Z , W is mod(Y,X), W \= 0, Q is X + 1,
	primechecker(Q,Y,Z).

isprime(2).
isprime(3).
isprime(P):- integer(P) , P > 3, P mod 2 =\= 0, \+ has_factor(P,3).
has_factor(N,L):- N mod L =:= 0.
has_factor(N,L):- L * L < N, L2 is L + 2, has_factor(N,L2).

prime_factors(N,L):- N > 0, prime_factors(N,L,2).
prime_factors(1,[],_):-!.
prime_factors(N,[F|L],F):- R is N // F, N =:= R * F, !, prime_factors(R,L,F).
prime_factors(N,L,F):- next_factor(N,F,NF),prime_factors(N,L,NF).
next_factor(_,2,3):-!.
next_factor(N,F,NF):- F * F < N, !, NF is F + 2.
next_factor(N,_,N).

prime_factors2(N,L):-prime_factors(N,X),theirEncode(X,L).

theirprime2(N,L):- N > 0, theirprime2(N,L,2).
theirprime2(1,[],_):-!.
theirprime2(N,[[F,M]|L],F):-divide(N,F,M,R),!,next_factor(R,F,NF),
	theirprime2(R,L,NF).
theirprime2(N,L,F):- !, next_factor(N,F,NF), theirprime2(N,L,NF).
divide(N,F,M,R):- divi(N,F,M,R,0), M > 0.
divi(N,F,M,R,K):- S is N // F, N =:= S * F, !,
	K1 is K + 1, divi(S,F,M,R,K1).
divi(N,_,M,N,M).

iseven(X):- X mod 2 =:= 0.

prime(X,Y,Z):- iseven(X),X1 is X - 1, primelist(X1,Y,Z).
prime(X,Y,Z):- primelist(X,Y,Z).
primelist(X,Y,[]):- Y =< X.
primelist(X,Y,[X|Z]):- isprime(X), X1 is X + 2, primelist(X1,Y,Z).
primelist(X,Y,Z):- X1 is X + 2, primelist(X1,Y,Z).

goldbach(X,Q):- Z is X - 3,  goldbach(X,3,Z,Q).
goldbach(X,Y,Z,[Y,Z]):- isprime(Y), isprime(Z), X =:= Y + Z.
goldbach(X,Y,Z,Q):- Y1 is Y + 2, Z1 is Z - 2,goldbach(X,Y1,Z1,Q).

goldbach_list(A,B):- goldbach_list(A,B,2).
goldbach_list(A,B,L):- A =< 4, ! ,g_list(4,B,L).
goldbach_list(A,B,L):- A1 is ((A+1) // 2) * 2, g_list(A1,B,L).

g_list(A,B,_):- A > B, !.
g_list(A,B,L):- goldbach(A,[P,Q]),print_goldbach(A,P,Q,L),A2 is A + 2, g_list(A2,B,L).

print_goldbach(A,P,Q,L):- P >= L, !, writef('%t = %t + %t' , [A,P,Q]), nl.
print_goldbach(_,_,_,_).

gcd(X,0,X):- X > 0.
gcd(X,Y,G):- Y > 0, Z is X mod Y, gcd(Y,Z,G).
:- arithmetic_function(gcd/2).

coprime(X,Y):- gcd(X,Y,1).

totient_phi(1,1) :- !.
totient_phi(M,Phi):- t_phi(M,Phi,1,0).

t_phi(M,Phi,M,Phi):-!.
t_phi(M,Phi,K,C):- K < M , coprime(K,M), !, C1 is C + 1, K1 is K + 1, t_phi(M,Phi,K1,C1).
t_phi(M,Phi,K,C):- K < M , K1 is K + 1, t_phi(M,Phi,K1,C).

totient_phi2(M,Phi):-  prime_factors2(M,Z), listplace(Z,X), multlist(X,Phi).



listplace([],[]).
listplace([X|Tx],[Y|Ty]):-multnode(X,Y),listplace(Tx,Ty).
multnode([TP|P],Y):- P1 is P - 1, TP1 is TP - 1, Y is (P1 * (P ** TP1)).
multlist([N], N).
multlist([N|Ns], M) :-
  multlist(Ns, T), M is N * T.

and(A,B):- A,B.
or(A,_):- A.
or(_,B):- B.
equ(A,B):- or( and(A,B), and(not(A),not(B))).
xor(A,B):- not(equ(A,B)).
nor(A,B):- not(or(A,B)).
nand(A,B):- not(and(A,B)).
impl(A,B):- or(not(A),B).

bind(true).
bind(false).



