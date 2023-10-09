% | MATERIAL | PERMITTIVITY |
% |   FR-4   |      4.4     |
% |   DPI    |      -       |
% |   PPO    |      2.7*    |
% |   PTFE   |      2.1     |
% |   CEM-1  |      4.2     |
% |   CEM-2  |      4.1*    |
% |   CEM-3  |      4.0     |

dielectric(4) = struct('name', [], 'er', []);
dielectric(1).name = 'N/A';
dielectric(1).er = 9.8;
dielectric(2).name = 'FR4';
dielectric(2).er = 4.4;
dielectric(3).name = 'PPO';
dielectric(3).er = 2.7;
dielectric(4).name = 'PTFE';
dielectric(4).er = 2.1;
