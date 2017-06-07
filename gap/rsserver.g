################################################################
# interface to the Riemann Sphere Server

#TODO:
# -- implement all commands: button, arc, etc. as high-level

RSS := rec(f := fail, # the file descriptor
           server := fail, # the daemon: rec(address,port,client-url)
           queue := [], # messages waiting to be delivered
           session := fail, # the default session id
           window := fail, # the default window id
           readhookslot := fail, # the slot in OnCharReadHookInFuncs
           callback := rec() # functions to call when button pressed
           );

OnCharReadHookActive := true; # we'll use read hooks

################################################################ low level IO
RSS.startdaemon := function()
    local server, s, i;
    if RSS.server<>fail then
        Error("Server seems already started");
    fi;
    CHECKEXEC@FR("node");
    server := IO_Popen3(EXEC@FR.node,[Filename(Directory(PackageInfo("img")[1].InstallationPath),"rsserver/rsserver.js")]);
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
    
    CHECKEXEC@FR("browser",["chrome"],["/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"],["firefox"]);
    
    status := IO_fork();
    if status=0 then
        IO_execv(EXEC@FR.browser,[url]);
        IO_exit(-1); # shouldn't happen
    fi;
    if status < 0 then
        Error("Couldn't start the client ",EXEC@FR.browser);
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
        if arg=[] then return; fi; # already open, but we accept
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
    local s;
    s := StringXMLElement(rec(name := "downdata", attributes := a, content := r))[1];
    Info(InfoIMG,2,"Sending XML string ",s);
    IO_Write(RSS.f, s);
end;

RSS.enqueue := function()
    local c, s;
    c := IO_ReadLine(RSS.f);
    for s in ParseTreeXMLString(c).content do
        if s.name="PCDATA" then # ignore
        elif s.name="error" then
            Error("RSS server error: ",s);
        elif s.name="updata" then
            if IsBound(s.attributes.status) and s.attributes.status="button-click" and IsBound(RSS.callback.(s.attributes.object)) then
                RSS.callback.(s.attributes.object)(s);
            else
                Add(RSS.queue,s);
            fi;
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
    return Remove(RSS.queue,1);
end;

RSS.flush := function()
    repeat
        RSS.queue := [];
        RSS.recv(false);
    until RSS.queue = [];
end;

RSS.ack := function(status, message)
    local i;
    
    while true do
        for i in [1..Length(RSS.queue)] do
            if RSS.queue[i].attributes.status=status then
                return Remove(RSS.queue,i);
            fi;
        od;
        # queue exhausted, wait.
        RSS.enqueue();
    od;
end;

################################################################ data handling
COMPLEX2XML@ := function(c)
    local a;
    a := rec(re := String(RealPart(c)), im := String(ImaginaryPart(c)));
    return rec(name := "cn", attributes := a, content := 0);
end;

P1POINT2XML@ := function(p)
    local a, c;
    if p=P1infinity then
        return rec(name := "cn", attributes := rec(name := "infinity"), content := 0);
    else
        return COMPLEX2XML@(P1Coordinate(p));
    fi;
end;

P1MAP2XML@ := function(f)
    local c, z, cycles;
    c := CoefficientsOfP1Map(f);
    z := PCDATAATTRACTINGCYCLES@IMG(POSTCRITICALPOINTS@IMG(f));
    cycles := Cycles(PermList(z{[1..Length(z)]}[2]+1),[1..Length(z)]);
    z := z{[1..Length(z)]}[1];
    
    #@ there was a typo, "nom" instead of "numer", which crashed node.
    return rec(name := "function",
               attributes := rec(degree := String(DegreeOfP1Map(f))),
               content := Concatenation([rec(name := "numer", attributes := rec(), content := List(c[1],COMPLEX2XML@)),
                       rec(name := "denom", attributes := rec(), content := List(c[2],COMPLEX2XML@))],
                       List(cycles,c->rec(name := "cycle", attributes := rec(), content := List(c,i->P1POINT2XML@(z[i]))))));
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

RSS.newobject := function(type,arg...)
    local status, attr;
    if arg=[] then attr := rec(); else attr := arg[1]; fi;
    RSS.send(rec(session:=RSS.session,object:=RSS.window),[rec(name:=type,attributes:=attr,content:=[])]);
    status := RSS.ack("created","newobject");
    return status.content[1].attributes.id;    
end;

RSS.removeobject := function(oid)
    RSS.send(rec(session:=RSS.session,object:=oid,action:="remove"),[]);
    RSS.ack("removed","removeobject");
end;

RSS.populateobject := function(oid,content)
    RSS.send(rec(session:=RSS.session,object:=oid,action:="populate"),content);
    RSS.ack("updated","populateobject");
end;

################################################################ top level
RSS.newcanvas := function()
    local status, killbutton, canvas;
    killbutton := RSS.newobject("button",rec(name:="Close canvas"));
    canvas := RSS.newobject("canvas");
    RSS.callback.(killbutton) := function(s)
        RSS.removeobject(killbutton);
        RSS.removeobject(canvas);
    end;
    return canvas;
end;

RSS.putmap := function(cid,map)
    RSS.populateobject(cid,[P1MAP2XML@(map)]);
end;

RSS.putarc := function(cid,arc)
    RSS.populatoobject(cid,[]);
end;

RSS.putpoint := function(cid,point,arg...)
    local a, attr, content;
    attr := rec();
    content := [P1POINT2XML@(point)];
    for a in arg do
        if IsFloat(a) then
            attr.radius := a;
        elif IsString(a) then
            Add(content,rec(name:="label",attributes:=rec(),content:=a));
        fi;
    od;
    RSS.populateobject(cid,[rec(name:="point",attributes:=attr,content:=content)]);
end;

################################################################ tests
RSS.samplewindow := function(map)
    local cid;

    RSS.open();
    if IsP1Map(map) then
        map := P1MAP2XML@(map);
    elif IsInt(map) then
        map := rec(name:="function",attributes:=rec(type:="newton",degree:=String(map)),content:=[]);
    fi;
    cid := RSS.newobject(RSS.window,"canvas");
    RSS.populateobject(cid,[map]);
    RSS.putpoint(cid,P1Point(1.0,0.5));
end;
