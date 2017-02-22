################################################################
# interface to the Riemann Sphere Server

#TODO:
# -- acknowledge should be different from "create" message;
#    have RSS.ack(filters) do that
# -- make systematic use of RSS.session, sometimes of RSS.window
# -- make better use of hooks: they should also trigger functions
#    in case a button is pressed, e.g.
# -- implement all commands: button, arc, etc. as high-level

RSS := rec(f := fail, # the file descriptor
           server := fail, # the daemon: rec(address,port,client-url)
           queue := [], # messages waiting to be delivered
           session := fail, # the default session id
           window := fail, # the default window id
           readhookslot := fail # the slot in OnCharReadHookInFuncs
           );

OnCharReadHookActive := true; # we'll use read hooks

################################################################ low level IO
RSS.startdaemon := function()
    local server, s, i;
    if RSS.server<>fail then
        Error("Server seems already started");
    fi;
    server := IO_Popen3(IO_FindExecutable("node"),[Filename(Directory(PackageInfo("img")[1].InstallationPath),"rsserver/rsserver.js")]);
    for i in [1..3] do
        s := IO_ReadLine(server.stdout);
        if s<>"" then Remove(s); fi; # remove \n
        Info(InfoIMG,1,s);
        s := SplitString(s," ,");
        if i=2 then
            if Length(s)<5 or s{[1..5]}<>["WebSocket","server","is","running.","Type"] then
                Error("Bad reply from daemon, no Websocket server");
            fi;
            server.client := s[6];
        elif i=3 then
            if Length(s)<5 or s{[1..5]}<>["TCP","server","is","running","at"] then
                Error("Bad reply from daemon, no TCP server");
            fi;
            server.address := s[6];
            server.port := Int(s[9]);
        fi;
    od;
    RSS.server := server;
end;

RSS.startclient := function(arg)
    local status, url, cmd;
    if arg=[] then
        if IsRecord(RSS.server) then
            url := RSS.server.client;
        else
            Error("Either start the daemon with RSS.startdaemon() or give a URL for the client-side of the daemon");
        fi;
    elif Length(arg)=1 and IsString(arg[1]) then
        url := arg[1];
    else
        Error("Use: RSS.startclient([url::string])");
    fi;
    
    #@ do something smarter to find the executable
    if POSITION_SUBSTRING(GAPInfo.Architecture,"darwin",0)=fail then
        cmd := "chrome";
    else
        cmd := "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome";
    fi;
    
    status := IO_fork();
    if status=0 then
        IO_execv(cmd,[url]);
        IO_exit(-1); # shouldn't happen
    fi;
    if status < 0 then
        Error("Couldn't start the client ",cmd);
    fi;
    
    status := RSS.recv(true);
    if status.attributes.status<>"created" then
        Error("Couldn't create browser session");
    fi;
    
    RSS.session := status.attributes.session;
    RSS.window := status.content[1].attributes.id;
end;

RSS.open := function(arg)
    local t, s, address, port;
    
    if arg=[] then
        if RSS.server=fail then
            RSS.startdaemon();
        fi;
        address := RSS.server.address;
        port := RSS.server.port;
    elif Length(arg)=2 and IsString(arg[1]) and IsInt(arg[2]) then
        address := arg[1];
        port := arg[2];
    else
        Error("Use: RSS.open([<address::string> <port::int>])");
    fi;

    if RSS.f <> fail then
        Error("Port seems already open");
    fi;
    
    t := IO_socket(IO.PF_INET,IO.SOCK_STREAM,"tcp");

    if IO_connect(t,IO_MakeIPAddressPort(address,port)) = fail then
        IO_close(t);
        Error("Could not connect to daemon");
        return fail;
    fi;
    
    RSS.f := IO_WrapFD(t,16384,false); # need buffered IO for IO_ReadLine
    RSS.queue := []; # messages waiting to be delivered
    
    RSS.readhookslot := Length(OnCharReadHookInFuncs)+1;
    Add(OnCharReadHookInFuncs, function(slotindex) RSS.enqueue(); end);
    Add(OnCharReadHookInFds, t);

    s := RSS.recv(true);
    
    if RSS.session=fail or s.content=0 then
        RSS.startclient();
    else
        s := s.content[1];
        RSS.session := s.attributes.id;
        RSS.window := s.content[1].attributes.id;
    fi;
end;

RSS.close := function()
    if not IsFile(RSS.f) then
        Error("Port seems already closed");
    fi;
    
    IO_Close(RSS.f);
    Remove(OnCharReadHookInFuncs, RSS.readhookslot);
    Remove(OnCharReadHookInFds, RSS.readhookslot);
    RSS.f := fail;
end;

RSS.send := function(a,r)
    IO_Write(RSS.f, StringXMLElement(rec(name := "downdata",
            attributes := a, content := r))[1]);
end;

RSS.enqueue := function()
    local c, s;
    c := IO_ReadLine(RSS.f);
    for s in ParseTreeXMLString(c).content do
        if s.name<>"PCDATA" then
            Add(RSS.queue,s);
        fi;
    od;
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
        RSS.enqueue();
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

RSS.samplewindow := function(map)
    local cid;
    if IsP1Map(map) then map := function2xml(map); fi;
    RSS.releasebutton(RSS.session,RSS.window);
    cid := RSS.newobject(RSS.session,RSS.window,"canvas");
    RSS.populateobject(RSS.session,cid,[map]);
    RSS.populateobject(RSS.session,cid,[rec(name:="point",attributes:=rec(),content:=[p1point2xml(P1Point(1.0,0.5))])]);
end;

newton := n->rec(name:="function",attributes:=rec(type:="newton",degree:=String(n)),content:=[]);

#inversebasilica := function2xml(CompositionP1Map(P1z^-1,P1z^2-1,P1z^-1));

# sample:
# RSS.open();
# RSS.samplewindow(P1z^2-1);
