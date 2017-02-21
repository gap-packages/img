RSS := rec(f := fail, server := fail);

################################################################ low level IO
RSS.serve := function()
    local s, i;
    if RSS.server<>fail then
        Error("Server seems already started");
    fi;
    RSS.server := IO_Popen3(IO_FindExecutable("node"),[Filename(Directory(PackageInfo("img")[1].InstallationPath),"rsserver/rsserver.js")]);
    for i in [1..3] do
        s := SplitString(IO_ReadLine(RSS.server.stdout),"\n")[1];
        Info(InfoIMG,1,s);
        s := SplitString(s," ,");
        if i=2 then
            if s{[1..5]}<>["WebSocket","server","is","running.","Type"] then
                Error("Bad reply from server");
            fi;
            RSS.server.client := s[6];
        elif i=3 then
            if s{[1..5]}<>["TCP","server","is","running","at"] then
                Error("Bad reply from server");
            fi;
            RSS.server.address := s[6];
            RSS.server.port := Int(s[9]);
        fi;
    od;
end;

RSS.client := function()
    local status;
    if IO_fork()=0 then
        #@ obviously the executable place should change
        IO_execv("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",["http://127.0.0.1:1729"]);
        IO_exit(-1);
    fi;
    status := RSS.recv(true);
    if status.attributes.status<>"created" then
        Error("Couldn't create browser session");
    fi;
    return rec(session:=status.attributes.session,window:=status.content[1].attributes.id);
end;

RSS.open := function(arg)
    local t, s, address, port;
    
    if arg=[] then
        if RSS.server=fail then
            RSS.serve();
        fi;
        address := RSS.server.address;
        port := RSS.server.port;
    elif Length(arg)=2 and IsString(arg[1]) and IsInt(arg[2]) then
        address := arg[1];
        port := arg[2];
    else
        Error("Use: RSS.open [<address> <port>]");
    fi;

    if RSS.f <> fail then
        Error("Port seems already open");
    fi;
    
    t := IO_socket(IO.PF_INET,IO.SOCK_STREAM,"tcp");

    if IO_connect(t,IO_MakeIPAddressPort(address,port)) = fail then
        IO_close(t);
        return fail;
    fi;
    
    RSS.f := IO_WrapFD(t,16384,false); # need buffered IO for IO_ReadLine
    RSS.queue := []; # messages waiting to be delivered
    
    s := RSS.recv(true);
    
    #@ what should we do in case multiple sessions are open?
    #@ should we store the session / window id in RSS so they can be used by default by the commands?
    if s.content=0 then
        return rec(); # no open session
    else
        s := s.content[1];
        return rec(session:=s.attributes.id,window:=s.content[1].attributes.id);
    fi;
end;

RSS.close := function()
    if not IsFile(RSS.f) then
        Error("Port seems already closed");
    fi;
    
    IO_Close(RSS.f);
    
    RSS.f := fail;
end;

RSS.send := function(a,r)
    IO_Write(RSS.f, StringXMLElement(rec(name := "downdata",
            attributes := a, content := r))[1]);
end;

RSS.recv := function(arg)
    local c, s, blocking;
    
    if arg=[] then
        blocking := false;
    elif Length(arg)=1 and IsBool(arg[1]) then
        blocking := arg[1];
    else
        Error("Use: RSS.recv [<blocking>]");
    fi;
        
    while RSS.queue=[] do
        if not (blocking or IO_HasData(RSS.f)) then return fail; fi;
        c := IO_ReadLine(RSS.f);
        for s in ParseTreeXMLString(c).content do
            if s.name<>"PCDATA" then
                Add(RSS.queue,s);
            fi;
        od;
    od;
    s := Remove(RSS.queue,1);
    if s.name="error" then
        Error("RSS server error: ",s);
    elif s.name="updata" then
        return s;
    else
        Error("RSS server unknown message: ",s);
    fi;
end;

################################################################ data handling
complex2xml := function(c)
    local a;
    a := rec(re := String(RealPart(c)), im := String(ImaginaryPart(c)));
    return rec(name := "cn", attributes := a, content := 0);
end;

p1point2xml := function(p)
    local a, c;
    if p=P1infinity then
        return rec(name := "cn", attributes := rec(name := "infinity"), content := 0);
    else
        return complex2xml(P1Coordinate(p));
    fi;
end;

function2xml := function(f)
    local c, z, cycles;
    c := CoefficientsOfP1Map(f);
    z := PCDATAATTRACTINGCYCLES@IMG(POSTCRITICALPOINTS@IMG(f));
    cycles := Cycles(PermList(z{[1..Length(z)]}[2]+1),[1..Length(z)]);
    z := z{[1..Length(z)]}[1];
    
    #@ there was a typo, "nom" instead of "numer", which crashed node.
    return rec(name := "function",
               attributes := rec(degree := String(DegreeOfP1Map(f))),
               content := Concatenation([rec(name := "numer", attributes := rec(), content := List(c[1],complex2xml)),
                       rec(name := "denom", attributes := rec(), content := List(c[2],complex2xml))],
                       List(cycles,c->rec(name := "cycle", attributes := rec(), content := List(c,i->p1point2xml(z[i]))))));
end;

################################################################ interface
RSS.getsession := function(arg)
    local cmd;
    cmd := rec(action:="request");
    if Length(arg)>2 or not ForAll(arg,IsString) then
        Error("Use: RSS.request [<session-id> [<object-id>] ]");
    fi;
    if Length(arg)>=1 then
        cmd.session := arg[1];
    fi;
    if Length(arg)>=2 then
        cmd.object := arg[2];
    fi;
    RSS.send(cmd,[]);
    return RSS.recv(true);
end;

RSS.newobject := function(sid,oid,type)
    local status;
    RSS.send(rec(session:=sid,object:=oid),[rec(name:=type,attributes:=rec(),content:=[])]);
    status := RSS.recv(true);
    if status.attributes.status<>"created" then
        Error("Couldn't create new object");
    fi;
    return status.content[1].attributes.id;    
end;

RSS.removeobject := function(sid,oid)
    local status;
    RSS.send(rec(session:=sid,object:=oid,action:="remove"),[]);
    status := RSS.recv(true);
    if status.attributes.status<>"removed" then
        Error("Couldn't remove object");
    fi;
end;

RSS.populateobject := function(sid,oid,content)
    local status;
    RSS.send(rec(session:=sid,object:=oid,action:="populate"),content);
    status := RSS.recv(true);
    if status.attributes.status<>"updated" then
        Error("Couldn't populate object");
    fi;
end;

################################################################ some tests
RSS.releasebutton := function(sid,oid)
    local status;
    RSS.send(rec(session:=sid,object:=oid),[rec(name:="button",attributes:=rec(name:="Release ReadLine"),content:=[])]);
    status := RSS.recv(true);
    if status.attributes.status<>"created" then
        Error("Couldn't create button");
    fi;
    return status.content[1].attributes.id;
end;

RSS.samplewindow := function(sid,wid,map)
    local cid;
    RSS.releasebutton(sid,wid);
    cid := RSS.newobject(sid,wid,"canvas");
    RSS.populateobject(sid,cid,[map]);
    RSS.populateobject(sid,cid,[rec(name:="point",attributes:=rec(),content:=[p1point2xml(P1Point(1.0,0.5))])]);
end;

#RSS.serve();
#RSS.open("localhost",1728);

newton := n->rec(name:="function",attributes:=rec(type:="newton",degree:=String(n)),content:=[]);

basilica := function2xml(P1z^2-1);

inversebasilica := function2xml(CompositionP1Map(P1z^-1,P1z^2-1,P1z^-1));

basilica1234 := fail;

func13 := fail;

# sample:
# RSS.serve();
# RSS.
