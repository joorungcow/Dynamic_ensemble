function kb = init_keyboard(kb)
% Setup universal Mac/PC keyboard and keynames to OSX naming scheme
KbName('UnifyKeyNames');

kb.escKey = KbName('ESCAPE');
kb.oneKey = [KbName('1!'), KbName('1')];
kb.twoKey = [KbName('2@'), KbName('2')];
kb.threeKey = [KbName('3#'), KbName('3')];
kb.fourKey = [KbName('4$'), KbName('4')];
kb.fiveKey = [KbName('5%'), KbName('5')];
kb.sixKey = [KbName('6^'), KbName('6')];
kb.sevenKey = [KbName('7&'), KbName('7')];
kb.eightKey = [KbName('8*'), KbName('8')];
kb.nineKey = [KbName('9('), KbName('9')];
kb.zeroKey = [KbName('0)'), KbName('0')];

kb.upKey = KbName('UpArrow');
kb.downKey = KbName('DownArrow');
kb.leftKey = KbName('LeftArrow');
kb.rightKey = KbName('rightArrow');
kb.spaceKey = KbName('space');

kb.qKey = KbName('q');
kb.wKey = KbName('w');
kb.eKey = KbName('e');
kb.rKey = KbName('r');
kb.pKey = KbName('p');
kb.oKey = KbName('o');
kb.iKey = KbName('i');
kb.kKey = KbName('k');
kb.lKey = KbName('l');
kb.zKey = KbName('z');
kb.xKey = KbName('x');
kb.cKey = KbName('c');
kb.aKey = KbName('a');
kb.sKey = KbName('s');
kb.dKey = KbName('d');
kb.fKey = KbName('f');
kb.nKey = KbName('n');

% Buttonbox keys
kb.tKey = KbName('t');
kb.yKey = KbName('y');
kb.bKey = KbName('b');
kb.gKey = KbName('g');
kb.rKey = KbName('r');

if ~exist('kb','var')
    kb.useKbQueueCheck = 0;
end

if ~isfield(kb,'useKbQueueCheck')
    kb.useKbQueueCheck = 0;
end

% Initialize Keyboards
if ~kb.useKbQueueCheck
    [kb.keyIsDown, kb.secs, kb.keyCode] = KbCheck(-1);
    oldenablekeys = RestrictKeysForKbCheck([kb.escKey, kb.leftKey, kb.rightKey, kb.spaceKey]);
else
    kb.devices = PsychHID('Devices');
    kb.keysOfInterest=zeros(1,256);
%     kb.keysOfInterest(KbName({'r', 'g', 'b', 'y', 't', 'a'})) = 1;
    kb.keysOfInterest(KbName(kb.keyList)) = 1;
    KbQueueCreate([], kb.keysOfInterest);
    % Perform some other initializations
    KbQueueStart;
    KbQueueCheck;
end

% Initialize Sound
InitializePsychSound;
Snd('Close');